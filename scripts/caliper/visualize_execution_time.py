#!/usr/bin/env python3
"""
visualize_execution_time.py
Visualizes function execution times from Caliper output files
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

def plot_execution_time(df, top_n=20, output_file='execution_time.png'):
    """Create horizontal bar chart of execution times"""
    df['function'] = df['HARD-Action'].fillna('') + '/' + df['HARD-Task'].fillna('')
    df['function'] = df['function'].str.strip('/').replace('', 'root')
    df['time_ms'] = df['time-ns'] / 1e6
    
    df_sorted = df.nlargest(top_n, 'time_ms')
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(df_sorted)))
    bars = ax.barh(df_sorted['function'], df_sorted['time_ms'], color=colors)
    
    ax.set_xlabel('Execution Time (ms)', fontsize=12)
    ax.set_ylabel('Function', fontsize=12)
    ax.set_title(f'Top {top_n} Functions by Execution Time', fontsize=14, fontweight='bold')
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3, axis='x')
    
    for i, (idx, row) in enumerate(df_sorted.iterrows()):
        ax.text(row['time_ms'] + row['time_ms']*0.01, i, 
                f"{row['time_ms']:.1f} ms ({int(row['count'])} calls)", 
                va='center', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Chart saved to {output_file}")
    
    return df_sorted

def plot_task_breakdown(df, top_n=20, output_file='task_breakdown.png'):
    """Create horizontal bar chart showing time distribution by task"""
    df_tasks = df.groupby('HARD-Task').agg({
        'time-ns': 'sum',
        'count': 'sum'
    }).reset_index()
    
    df_tasks = df_tasks[df_tasks['HARD-Task'].notna()]
    df_tasks['time_ms'] = df_tasks['time-ns'] / 1e6
    df_tasks = df_tasks.nlargest(top_n, 'time_ms')
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    colors = plt.cm.plasma(np.linspace(0, 1, len(df_tasks)))
    bars = ax.barh(df_tasks['HARD-Task'], df_tasks['time_ms'], color=colors)
    
    ax.set_xlabel('Execution Time (ms)', fontsize=12)
    ax.set_ylabel('Task', fontsize=12)
    ax.set_title(f'Top {top_n} Tasks by Execution Time', fontsize=14, fontweight='bold')
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3, axis='x')
    
    for i, (idx, row) in enumerate(df_tasks.iterrows()):
        ax.text(row['time_ms'] + row['time_ms']*0.01, i, 
                f"{row['time_ms']:.1f} ms ({int(row['count'])} calls)", 
                va='center', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Chart saved to {output_file}")

def plot_action_breakdown(df, top_n=20, output_file='action_breakdown.png'):
    """Create horizontal bar chart showing time distribution by action"""
    df_actions = df.groupby('HARD-Action').agg({
        'time-ns': 'sum',
        'count': 'sum'
    }).reset_index()
    
    df_actions = df_actions[df_actions['HARD-Action'].notna()]
    df_actions['time_ms'] = df_actions['time-ns'] / 1e6
    df_actions = df_actions.nlargest(top_n, 'time_ms')
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    colors = plt.cm.coolwarm(np.linspace(0, 1, len(df_actions)))
    bars = ax.barh(df_actions['HARD-Action'], df_actions['time_ms'], color=colors)
    
    ax.set_xlabel('Execution Time (ms)', fontsize=12)
    ax.set_ylabel('Action', fontsize=12)
    ax.set_title(f'Top {top_n} Actions by Execution Time', fontsize=14, fontweight='bold')
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3, axis='x')
    
    for i, (idx, row) in enumerate(df_actions.iterrows()):
        ax.text(row['time_ms'] + row['time_ms']*0.01, i, 
                f"{row['time_ms']:.1f} ms ({int(row['count'])} calls)", 
                va='center', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Chart saved to {output_file}")

def print_statistics(df):
    """Print execution time statistics"""
    print("\n" + "="*70)
    print("EXECUTION TIME STATISTICS")
    print("="*70)
    
    total_time_ms = df['time-ns'].sum() / 1e6
    total_time_s = total_time_ms / 1000
    print(f"\nTotal execution time: {total_time_ms:.2f} ms ({total_time_s:.3f} s)")
    print(f"Total function calls: {int(df['count'].sum())}")
    print(f"Number of functions: {len(df)}")
    
    df['function'] = df['HARD-Action'].fillna('') + '/' + df['HARD-Task'].fillna('')
    df['function'] = df['function'].str.strip('/').replace('', 'root')
    df['time_ms'] = df['time-ns'] / 1e6
    
    print("\n--- Top 10 Most Expensive Functions ---")
    top10 = df.nlargest(10, 'time_ms')
    print(f"{'Function':<50} {'Time (ms)':>12} {'% Total':>10} {'Calls':>10}")
    print("-" * 85)
    
    for idx, row in top10.iterrows():
        func_name = row['function'][:47] + '...' if len(row['function']) > 50 else row['function']
        print(f"{func_name:<50} {row['time_ms']:>12.2f} {row['time_ms']/total_time_ms*100:>9.1f}% {int(row['count']):>10}")

def main():
    parser = argparse.ArgumentParser(
        description='Visualize execution times from Caliper output files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python visualize_execution_time.py data.cali
  python visualize_execution_time.py data.cali -n 30
  python visualize_execution_time.py data.cali --top-n 15
        """
    )
    
    parser.add_argument('input_file', 
                        help='Caliper output file to analyze')
    parser.add_argument('-n', '--top-n', 
                        type=int, 
                        default=20,
                        help='Number of top functions to display (default: 20)')
    parser.add_argument('--output-prefix',
                        default='',
                        help='Prefix for output files')
    
    args = parser.parse_args()
    
    try:
        df = parse_caliper_file(args.input_file)
        print(f"Loaded {len(df)} entries from {args.input_file}")
        
        prefix = args.output_prefix + '_' if args.output_prefix else ''
        
        print_statistics(df)
        plot_execution_time(df, top_n=args.top_n, 
                          output_file=f'{prefix}execution_time.png')
        plot_task_breakdown(df, top_n=args.top_n, 
                          output_file=f'{prefix}task_breakdown.png')
        plot_action_breakdown(df, top_n=args.top_n, 
                            output_file=f'{prefix}action_breakdown.png')
        
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