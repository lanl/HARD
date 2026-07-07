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

def plot_execution_time(df, top_n=20, output_file='execution_time.png'):
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
            f"{row['time_ms']:.1f} ms ({int(row['count'])} calls)", va='center', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Chart saved to {output_file}")
    
    return df_sorted

def plot_task_breakdown(df, top_n=20, output_file='task_breakdown.png'):
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

    parser.add_argument('input_file', help='Caliper output file to analyze')
    parser.add_argument('-n', '--top-n', type=int, default=20,
                        help='Number of top functions to display (default: 20)')
    parser.add_argument('--output-prefix', default='',
                        help='Prefix for output files')
    
    args = parser.parse_args()
    
    df = parse_caliper_file(args.input_file)
    print(f"Loaded {len(df)} entries from {args.input_file}")
    
    prefix = args.output_prefix + '_' if args.output_prefix else ''
    
    plot_execution_time(df, top_n=args.top_n, output_file=f'{prefix}execution_time.png')
    plot_task_breakdown(df, top_n=args.top_n, output_file=f'{prefix}task_breakdown.png')
    plot_action_breakdown(df, top_n=args.top_n, output_file=f'{prefix}action_breakdown.png')
    
    print("\nAnalysis completed successfully!")

if __name__ == "__main__":
    main()