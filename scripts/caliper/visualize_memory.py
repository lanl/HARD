#!/usr/bin/env python3
"""
visualize_memory.py
Visualizes memory consumption from Caliper output files
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

def plot_memory_by_id(df, output_file='memory_by_id.png'):
    """Plot memory consumption by function ID (starting from 1)"""
    df['rss_mb'] = df['avg-rss-kb'] / 1024
    df['vm_mb'] = df['avg-vm-kb'] / 1024
    df['data_mb'] = df['avg-data-kb'] / 1024
    
    df_filtered = df[df.index > 0].copy()
    
    fig, ax = plt.subplots(figsize=(16, 8))
    
    x = df_filtered.index
    ax.plot(x, df_filtered['rss_mb'], 'b-', linewidth=2, label='RSS (Resident Set Size)', marker='o', markersize=3)
    ax.plot(x, df_filtered['vm_mb'], 'r-', linewidth=2, label='VM (Virtual Memory)', marker='s', markersize=3)
    ax.plot(x, df_filtered['data_mb'], 'g-', linewidth=2, label='Data Segment', marker='^', markersize=3)
    
    ax.fill_between(x, df_filtered['rss_mb'], alpha=0.3, color='blue')
    
    ax.set_xlabel('Function ID', fontsize=12)
    ax.set_ylabel('Memory (MB)', fontsize=12)
    ax.set_title('Memory Consumption Evolution', fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Chart saved to {output_file}")

def plot_task_memory(df, top_n=20, output_file='task_memory.png'):
    """Plot average memory consumption by task"""
    df_tasks = df.groupby('HARD-Task').agg({
        'avg-rss-kb': 'mean',
        'count': 'sum'
    }).reset_index()
    
    df_tasks = df_tasks[df_tasks['HARD-Task'].notna()]
    
    # Filter to keep only valid task names (containing letters/underscores, not just numbers)
    df_tasks = df_tasks[df_tasks['HARD-Task'].astype(str).str.contains(r'[a-zA-Z_]', regex=True)]
    
    df_tasks['avg_rss_mb'] = df_tasks['avg-rss-kb'] / 1024
    df_tasks = df_tasks.nlargest(top_n, 'avg_rss_mb')
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    colors = plt.cm.plasma(np.linspace(0, 1, len(df_tasks)))
    bars = ax.barh(df_tasks['HARD-Task'], df_tasks['avg_rss_mb'], color=colors)
    
    ax.set_xlabel('Average Memory RSS (MB)', fontsize=12)
    ax.set_ylabel('Task', fontsize=12)
    ax.set_title(f'Top {top_n} Tasks by Average Memory Consumption', fontsize=14, fontweight='bold')
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3, axis='x')
    
    for i, (idx, row) in enumerate(df_tasks.iterrows()):
        ax.text(row['avg_rss_mb'] + row['avg_rss_mb']*0.01, i, 
                f"{row['avg_rss_mb']:.1f} MB ({int(row['count'])} calls)", 
                va='center', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Chart saved to {output_file}")

def plot_time_vs_memory(df, output_file='time_vs_memory.png'):
    """Plot execution time vs memory consumption"""
    df['time_ms'] = df['time-ns'] / 1e6
    df['rss_mb'] = df['avg-rss-kb'] / 1024
    
    # Filter out entries with zero or very low execution time and valid memory
    df_filtered = df[(df['time_ms'] > 0.1) & (df['rss_mb'] > 0)].copy()
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    scatter = ax.scatter(df_filtered['time_ms'], df_filtered['rss_mb'], 
                        alpha=0.6, s=50, c='steelblue', edgecolors='black', linewidth=0.5)
    
    ax.set_xlabel('Execution Time (ms)', fontsize=12)
    ax.set_ylabel('Average Memory RSS (MB)', fontsize=12)
    ax.set_title('Execution Time vs Memory Consumption', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    if len(df_filtered) > 0:
        z = np.polyfit(df_filtered['time_ms'], df_filtered['rss_mb'], 1)
        p = np.poly1d(z)
        ax.plot(df_filtered['time_ms'].sort_values(), 
                p(df_filtered['time_ms'].sort_values()), 
                "r--", alpha=0.5, linewidth=2, label=f'Trend line (slope={z[0]:.2e})')
        ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Chart saved to {output_file}")
    
def print_memory_statistics(df):
    """Print memory consumption statistics"""
    print("\n" + "="*70)
    print("MEMORY CONSUMPTION STATISTICS")
    print("="*70)
    
    max_rss = df['max-rss-kb'].max() / 1024
    avg_rss = df['avg-rss-kb'].mean() / 1024
    max_vm = df['max-vm-kb'].max() / 1024
    
    print(f"\nMaximum RSS memory: {max_rss:.2f} MB")
    print(f"Global average RSS memory: {avg_rss:.2f} MB")
    print(f"Maximum VM memory: {max_vm:.2f} MB")
    
    df_tasks = df.groupby('HARD-Task').agg({
        'avg-rss-kb': 'mean'
    }).reset_index()
    df_tasks = df_tasks[df_tasks['HARD-Task'].notna()]
    df_tasks['avg_rss_mb'] = df_tasks['avg-rss-kb'] / 1024
    
    print("\n--- Top 10 Tasks by Average Memory Consumption ---")
    top10 = df_tasks.nlargest(10, 'avg_rss_mb')
    print(f"{'Task':<50} {'Avg RSS (MB)':>15}")
    print("-" * 70)
    
    for idx, row in top10.iterrows():
        task_name = row['HARD-Task'][:47] + '...' if len(str(row['HARD-Task'])) > 50 else row['HARD-Task']
        print(f"{task_name:<50} {row['avg_rss_mb']:>15.2f}")

def main():
    parser = argparse.ArgumentParser(
        description='Visualize memory consumption from Caliper output files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python visualize_memory.py data.cali
  python visualize_memory.py data.cali -n 30
  python visualize_memory.py data.cali --top-n 15 --output-prefix baseline
        """
    )
    
    parser.add_argument('input_file', 
                        help='Caliper output file to analyze')
    parser.add_argument('-n', '--top-n', 
                        type=int, 
                        default=20,
                        help='Number of top tasks to display (default: 20)')
    parser.add_argument('--output-prefix',
                        default='',
                        help='Prefix for output files')
    
    args = parser.parse_args()
    
    try:
        df = parse_caliper_file(args.input_file)
        print(f"Loaded {len(df)} entries from {args.input_file}")
        
        prefix = args.output_prefix + '_' if args.output_prefix else ''
        
        print_memory_statistics(df)
        plot_memory_by_id(df, output_file=f'{prefix}memory_by_id.png')
        plot_task_memory(df, top_n=args.top_n, output_file=f'{prefix}task_memory.png')
        plot_time_vs_memory(df, output_file=f'{prefix}time_vs_memory.png')
        
        print("\n✅ Analysis completed successfully!")
        
    except FileNotFoundError:
        print(f"Error: File '{args.input_file}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()