#!/usr/bin/env python3
"""
compare_caliper_outputs.py
Compare two Caliper output files to identify performance improvements/regressions
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse

import matplotlib
matplotlib.use("Agg")

def parse_caliper_file(filename):
    """Parse Caliper file and return DataFrame"""
    df = pd.read_csv(filename, delim_whitespace=True, skiprows=0)
    return df

def compare_files(file1, file2, label1="Version 1", label2="Version 2"):
    """Compare two Caliper files"""
    df1 = parse_caliper_file(file1)
    df2 = parse_caliper_file(file2)
    
    print("\n" + "="*70)
    print(f"COMPARISON: {label1} vs {label2}")
    print("="*70)
    
    total_time1 = df1['time-ns'].sum() / 1e9
    total_time2 = df2['time-ns'].sum() / 1e9
    
    print(f"\n{label1}:")
    print(f"  Total time: {total_time1:.3f} s")
    print(f"  Total calls: {int(df1['count'].sum())}")
    print(f"  Max RSS memory: {df1['max-rss-kb'].max() / 1024:.2f} MB")
    
    print(f"\n{label2}:")
    print(f"  Total time: {total_time2:.3f} s")
    print(f"  Total calls: {int(df2['count'].sum())}")
    print(f"  Max RSS memory: {df2['max-rss-kb'].max() / 1024:.2f} MB")
    
    time_diff = ((total_time2 - total_time1) / total_time1) * 100 if total_time1 > 0 else 0
    mem_diff = ((df2['max-rss-kb'].max() - df1['max-rss-kb'].max()) / df1['max-rss-kb'].max()) * 100 if df1['max-rss-kb'].max() > 0 else 0
    
    print(f"\n--- DIFFERENCES ---")
    print(f"Execution time: {time_diff:+.2f}% {'🚀 IMPROVEMENT' if time_diff < 0 else '⚠️  REGRESSION' if time_diff > 0 else '➖ IDENTICAL'}")
    print(f"Max memory: {mem_diff:+.2f}% {'✅ IMPROVEMENT' if mem_diff < 0 else '⚠️  REGRESSION' if mem_diff > 0 else '➖ IDENTICAL'}")
    
    comparison = compare_tasks(df1, df2, label1, label2)
    
    return df1, df2, comparison

def compare_tasks(df1, df2, label1, label2):
    """Compare tasks between two versions"""
    df1_tasks = df1.groupby('HARD-Task').agg({
        'time-ns': 'sum',
        'count': 'sum'
    }).reset_index()
    
    df2_tasks = df2.groupby('HARD-Task').agg({
        'time-ns': 'sum',
        'count': 'sum'
    }).reset_index()
    
    df1_tasks = df1_tasks[df1_tasks['HARD-Task'].notna()]
    df2_tasks = df2_tasks[df2_tasks['HARD-Task'].notna()]
    
    df1_tasks = df1_tasks[df1_tasks['HARD-Task'].astype(str).str.contains(r'[a-zA-Z_]', regex=True)]
    df2_tasks = df2_tasks[df2_tasks['HARD-Task'].astype(str).str.contains(r'[a-zA-Z_]', regex=True)]
    
    merged = pd.merge(df1_tasks, df2_tasks, on='HARD-Task', suffixes=('_v1', '_v2'), how='outer')
    
    merged['time_ms_v1'] = merged['time-ns_v1'] / 1e6
    merged['time_ms_v2'] = merged['time-ns_v2'] / 1e6
    
    merged['time_diff_pct'] = ((merged['time_ms_v2'] - merged['time_ms_v1']) / 
                               merged['time_ms_v1'].replace(0, np.nan)) * 100
    merged['time_diff_abs_pct'] = merged['time_diff_pct'].abs()
    
    print("\n--- TOP 10 IMPROVEMENTS (tasks) ---")
    improvements = merged[merged['time_diff_pct'] < -0.1].nsmallest(10, 'time_diff_pct')
    if len(improvements) == 0:
        print("No significant improvements detected")
    else:
        for idx, row in improvements.iterrows():
            print(f"{row['HARD-Task'][:50]:50s} {row['time_diff_pct']:+7.2f}% "
                  f"({row['time_ms_v1']:.2f} → {row['time_ms_v2']:.2f} ms)")
    
    print("\n--- TOP 10 REGRESSIONS (tasks) ---")
    regressions = merged[merged['time_diff_pct'] > 0.1].nlargest(10, 'time_diff_pct')
    if len(regressions) == 0:
        print("No significant regressions detected")
    else:
        for idx, row in regressions.iterrows():
            print(f"{row['HARD-Task'][:50]:50s} {row['time_diff_pct']:+7.2f}% "
                  f"({row['time_ms_v1']:.2f} → {row['time_ms_v2']:.2f} ms)")
    
    return merged

def plot_comparison(df1, df2, comparison, top_n=20, label1="Version 1", label2="Version 2", 
                   output_file='comparison.png'):
    """Create comparison charts"""
    fig = plt.figure(figsize=(18, 12))
    
    # Chart 1: Total execution time comparison
    ax1 = plt.subplot(2, 2, 1)
    total_times = [df1['time-ns'].sum() / 1e6, df2['time-ns'].sum() / 1e6]
    colors = ['red' if total_times[1] > total_times[0] else 'green', 'steelblue']
    bars = ax1.bar([label1, label2], total_times, color=colors, alpha=0.7)
    ax1.set_ylabel('Total Time (ms)', fontsize=12)
    ax1.set_title('Total Execution Time', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')
    
    for bar, time in zip(bars, total_times):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{time:.1f} ms',
                ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Chart 2: Maximum memory comparison
    ax2 = plt.subplot(2, 2, 2)
    max_mems = [df1['max-rss-kb'].max() / 1024, df2['max-rss-kb'].max() / 1024]
    colors = ['red' if max_mems[1] > max_mems[0] else 'green', 'steelblue']
    bars = ax2.bar([label1, label2], max_mems, color=colors, alpha=0.7)
    ax2.set_ylabel('Max Memory (MB)', fontsize=12)
    ax2.set_title('Maximum Memory Consumption', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    
    for bar, mem in zip(bars, max_mems):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{mem:.1f} MB',
                ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Chart 3: Top tasks comparison by execution time
    ax3 = plt.subplot(2, 2, 3)
    
    top_tasks = comparison.nlargest(top_n, 'time_ms_v1')
    x = np.arange(len(top_tasks))
    width = 0.35
    
    bars1 = ax3.barh(x + width/2, top_tasks['time_ms_v1'], width, label=label1, alpha=0.8, color='steelblue')
    bars2 = ax3.barh(x - width/2, top_tasks['time_ms_v2'], width, label=label2, alpha=0.8, color='orange')
    
    ax3.set_xlabel('Execution Time (ms)', fontsize=12)
    ax3.set_ylabel('Task', fontsize=12)
    ax3.set_title(f'Top {top_n} Tasks: Execution Time Comparison', fontsize=14, fontweight='bold')
    ax3.set_yticks(x)
    ax3.set_yticklabels(top_tasks['HARD-Task'], fontsize=8)
    ax3.legend(loc='best')
    ax3.grid(True, alpha=0.3, axis='x')
    ax3.invert_yaxis()
    
    # Chart 4: Tasks with largest absolute percentage change
    ax4 = plt.subplot(2, 2, 4)
    
    comp_filtered = comparison[comparison['time_diff_pct'].notna()].copy()
    top_changes = comp_filtered.nlargest(top_n, 'time_diff_abs_pct')
    
    colors_diff = ['green' if x < 0 else 'red' for x in top_changes['time_diff_pct']]
    bars = ax4.barh(range(len(top_changes)), top_changes['time_diff_pct'], color=colors_diff, alpha=0.7)
    ax4.set_yticks(range(len(top_changes)))
    ax4.set_yticklabels(top_changes['HARD-Task'], fontsize=8)
    ax4.set_xlabel('Performance Change (%)', fontsize=12)
    ax4.set_ylabel('Task', fontsize=12)
    ax4.set_title(f'Top {top_n} Tasks by Absolute Performance Change\n(green=improvement, red=regression)', 
                  fontsize=14, fontweight='bold')
    ax4.axvline(x=0, color='black', linestyle='--', linewidth=1)
    ax4.grid(True, alpha=0.3, axis='x')
    ax4.invert_yaxis()
    
    for i, (idx, row) in enumerate(top_changes.iterrows()):
        value = row['time_diff_pct']
        x_pos = value + (2 if value > 0 else -2)
        ax4.text(x_pos, i, f'{value:+.1f}%', 
                va='center', ha='left' if value > 0 else 'right', 
                fontsize=8, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nComparison chart saved to {output_file}")

def generate_report(df1, df2, comparison, label1, label2, output_file='comparison_report.txt'):
    """Generate detailed text report"""
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("="*70 + "\n")
        f.write(f"HPC PERFORMANCE COMPARISON REPORT - {label1} vs {label2}\n")
        f.write("="*70 + "\n\n")
        
        total_time1 = df1['time-ns'].sum() / 1e9
        total_time2 = df2['time-ns'].sum() / 1e9
        time_diff = ((total_time2 - total_time1) / total_time1) * 100 if total_time1 > 0 else 0
        speedup = total_time1 / total_time2 if total_time2 > 0 else 1
        
        mem1 = df1['max-rss-kb'].max() / 1024
        mem2 = df2['max-rss-kb'].max() / 1024
        mem_diff = ((mem2 - mem1) / mem1) * 100 if mem1 > 0 else 0
        
        f.write("EXECUTIVE SUMMARY\n")
        f.write("-" * 70 + "\n\n")
        
        f.write(f"Execution time:\n")
        f.write(f"  {label1}: {total_time1:.3f} s\n")
        f.write(f"  {label2}: {total_time2:.3f} s\n")
        f.write(f"  Difference: {time_diff:+.2f}%\n")
        f.write(f"  Speedup: {speedup:.3f}x\n\n")
        
        f.write(f"Memory consumption:\n")
        f.write(f"  {label1}: {mem1:.2f} MB\n")
        f.write(f"  {label2}: {mem2:.2f} MB\n")
        f.write(f"  Difference: {mem_diff:+.2f}%\n\n")
        
        f.write("VERDICT: ")
        if time_diff < -5 and mem_diff < 5:
            f.write("✅ SIGNIFICANT IMPROVEMENT\n")
        elif time_diff < 0:
            f.write("✅ Slight improvement\n")
        elif abs(time_diff) < 0.5:
            f.write("➖ Similar performance\n")
        else:
            f.write("⚠️  REGRESSION\n")
        
        f.write("\n" + "="*70 + "\n\n")
        f.write("TOP 20 IMPROVEMENTS (tasks)\n")
        f.write("-" * 70 + "\n")
        f.write(f"{'Task':<50} {'Diff %':>10} {'V1 (ms)':>12} {'V2 (ms)':>12}\n")
        f.write("-" * 70 + "\n")
        
        improvements = comparison[comparison['time_diff_pct'] < -0.1].nsmallest(20, 'time_diff_pct')
        if len(improvements) == 0:
            f.write("No significant improvements detected\n")
        else:
            for idx, row in improvements.iterrows():
                f.write(f"{row['HARD-Task'][:50]:<50} {row['time_diff_pct']:>9.2f}% "
                       f"{row['time_ms_v1']:>11.2f}  {row['time_ms_v2']:>11.2f}\n")
        
        f.write("\n" + "="*70 + "\n\n")
        f.write("TOP 20 REGRESSIONS (tasks)\n")
        f.write("-" * 70 + "\n")
        f.write(f"{'Task':<50} {'Diff %':>10} {'V1 (ms)':>12} {'V2 (ms)':>12}\n")
        f.write("-" * 70 + "\n")
        
        regressions = comparison[comparison['time_diff_pct'] > 0.1].nlargest(20, 'time_diff_pct')
        if len(regressions) == 0:
            f.write("No significant regressions detected\n")
        else:
            for idx, row in regressions.iterrows():
                f.write(f"{row['HARD-Task'][:50]:<50} {row['time_diff_pct']:>9.2f}% "
                       f"{row['time_ms_v1']:>11.2f}  {row['time_ms_v2']:>11.2f}\n")
        
        f.write("\n" + "="*70 + "\n\n")
        f.write("TOP 20 TASKS BY ABSOLUTE PERFORMANCE CHANGE\n")
        f.write("-" * 70 + "\n")
        f.write(f"{'Task':<50} {'Diff %':>10} {'V1 (ms)':>12} {'V2 (ms)':>12}\n")
        f.write("-" * 70 + "\n")
        
        top_changes = comparison[comparison['time_diff_pct'].notna()].nlargest(20, 'time_diff_abs_pct')
        for idx, row in top_changes.iterrows():
            f.write(f"{row['HARD-Task'][:50]:<50} {row['time_diff_pct']:>9.2f}% "
                   f"{row['time_ms_v1']:>11.2f}  {row['time_ms_v2']:>11.2f}\n")
        
        f.write("\n" + "="*70 + "\n")
        f.write("END OF REPORT\n")
        f.write("="*70 + "\n")
    
    print(f"Detailed report saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Compare two Caliper output files for performance analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python compare_caliper_outputs.py baseline.cali optimized.cali
  python compare_caliper_outputs.py file1.cali file2.cali -n 30
  python compare_caliper_outputs.py baseline.cali opt.cali --label1 "Baseline" --label2 "Optimized"
        """
    )
    
    parser.add_argument('file1', help='First Caliper output file (baseline)')
    parser.add_argument('file2', help='Second Caliper output file (comparison)')
    parser.add_argument('-n', '--top-n', type=int, default=20,
                        help='Number of top tasks to display (default: 20)')
    parser.add_argument('--label1', default='Version 1',
                        help='Label for first file')
    parser.add_argument('--label2', default='Version 2',
                        help='Label for second file')
    parser.add_argument('--output-prefix', default='comparison',
                        help='Prefix for output files')
    
    args = parser.parse_args()
    
    try:
        df1, df2, comparison = compare_files(args.file1, args.file2, args.label1, args.label2)
        
        plot_comparison(df1, df2, comparison, top_n=args.top_n, 
                       label1=args.label1, label2=args.label2,
                       output_file=f'{args.output_prefix}.png')
        
        generate_report(df1, df2, comparison, args.label1, args.label2,
                       output_file=f'{args.output_prefix}_report.txt')
        
        print("\n✅ Analysis completed successfully!")
        print(f"\nGenerated files:")
        print(f"  - {args.output_prefix}.png (charts)")
        print(f"  - {args.output_prefix}_report.txt (detailed report)")
        
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()