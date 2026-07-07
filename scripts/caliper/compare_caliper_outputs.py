#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse

import matplotlib
matplotlib.use("Agg")

def parse_caliper_file(filename):
    df = pd.read_csv(filename, delim_whitespace=True, skiprows=0)
    return df

def compare_files(file1, file2, label1="Version 1", label2="Version 2"):
    df1 = parse_caliper_file(file1)
    df2 = parse_caliper_file(file2)
    
    total_time1 = df1['time-ns'].sum() / 1e9
    total_time2 = df2['time-ns'].sum() / 1e9
    
    time_diff = ((total_time2 - total_time1) / total_time1) * 100 if total_time1 > 0 else 0
    mem_diff = ((df2['max-rss-kb'].max() - df1['max-rss-kb'].max()) / df1['max-rss-kb'].max()) * 100 if df1['max-rss-kb'].max() > 0 else 0
    
    comparison = compare_tasks(df1, df2, label1, label2)
    
    return df1, df2, comparison

def compare_tasks(df1, df2, label1, label2):
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
    
    improvements = merged[merged['time_diff_pct'] < -0.1].nsmallest(10, 'time_diff_pct')
    regressions = merged[merged['time_diff_pct'] > 0.1].nlargest(10, 'time_diff_pct')
    
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
    
    df1, df2, comparison = compare_files(args.file1, args.file2, args.label1, args.label2)
    
    plot_comparison(df1, df2, comparison, top_n=args.top_n, 
                    label1=args.label1, label2=args.label2,
                    output_file=f'{args.output_prefix}.png')
            
    print("\nAnalysis completed successfully!")
    print(f"\nGenerated files:")
    print(f"  - {args.output_prefix}.png (charts)")

if __name__ == "__main__":
    main()