#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

def parse_worker_count_from_arg(arg):
    if '=' in arg:
        parts = arg.rsplit('=', 1)
        filename = parts[0]
        try:
            worker_count = int(parts[1])
            return filename, worker_count
        except ValueError:
            print(f"Warning: Invalid worker count in '{arg}', skipping...")
            return None, None
    else:
        print(f"Warning: Argument '{arg}' missing '=X' format, skipping...")
        return None, None

def parse_caliper_file(filename):
    df = pd.read_csv(filename, delim_whitespace=True, skiprows=0)
    return df

def load_all_files(file_args):
    data = {}
    
    for arg in file_args:
        filename, worker_count = parse_worker_count_from_arg(arg)
        
        if filename is None or worker_count is None:
            continue
        
        try:
            df = parse_caliper_file(filename)
            data[worker_count] = df
            print(f"Loaded {filename}: {worker_count} workers, {len(df)} entries")
        except Exception as e:
            print(f"Error loading {filename}: {e}")
            continue
    
    if not data:
        raise ValueError("No valid files loaded!")
    
    return data

def get_task_scaling_data(data_dict):
    """Extract task-level scaling data using task ID as unique identifier"""
    all_tasks = {}
    task_time_percentages = {}
    task_id_to_name = {}
    
    for worker_count, df in sorted(data_dict.items()):
        df_filtered = df[df['HARD-Task'].notna()].copy()
        df_filtered = df_filtered[df_filtered['HARD-Task'].astype(str).str.contains(r'[a-zA-Z_]', regex=True)]
        df_filtered['time_ms'] = df_filtered['time-ns'] / 1e6
        
        for _, row in df_filtered.iterrows():
            task_id = row['id']
            task_name = row['HARD-Task']
            task_key = (task_id, task_name)
            
            # Stocker le mapping
            if task_id not in task_id_to_name:
                task_id_to_name[task_id] = task_name
            
            if task_key not in all_tasks:
                all_tasks[task_key] = {}
            all_tasks[task_key][worker_count] = row['time_ms']
            
            if worker_count == min(data_dict.keys()):
                task_time_percentages[task_key] = row['total-time-%']
    
    return all_tasks, task_time_percentages

def plot_task_scaling(all_tasks, data_dict, top_n=10, output_file='scaling_tasks.png'):
    # Get minimum worker count
    min_workers = min(data_dict.keys())
    
    # Find top N tasks for minimum worker count
    tasks_at_min = [(task_key, times.get(min_workers, 0)) 
                    for task_key, times in all_tasks.items() 
                    if min_workers in times]
    tasks_at_min.sort(key=lambda x: x[1], reverse=True)
    top_tasks = [task_key for task_key, _ in tasks_at_min[:top_n]]
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    colors = plt.cm.tab20(np.linspace(0, 1, len(top_tasks)))
    
    for i, task_key in enumerate(top_tasks):
        task_id, task_name = task_key
        task_data = all_tasks[task_key]
        workers = sorted(task_data.keys())
        times = [task_data[w] for w in workers]
        
        label = f"{task_name[:35]} (ID:{task_id})"
        ax.plot(workers, times, marker='o', linewidth=2, markersize=8,
                label=label, color=colors[i])
    
    ax.set_xlabel('Number of Workers', fontsize=12)
    ax.set_ylabel('Execution Time (ms)', fontsize=12)
    ax.set_title(f'Scaling Performance: Top {top_n} Tasks\n(Based on execution time with {min_workers} workers)',
                 fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')

def plot_total_scaling(data_dict, output_file='scaling_total.png'):
    workers = []
    total_times = []
    
    for worker_count, df in sorted(data_dict.items()):
        workers.append(worker_count)
        total_time_ms = df['time-ns'].sum() / 1e6
        total_times.append(total_time_ms)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Chart 1: Total execution time
    ax1.plot(workers, total_times, marker='o', linewidth=3, markersize=10,
             color='steelblue', label='Actual')
    
    baseline_time = total_times[0]
    baseline_workers = workers[0]
    ideal_times = [baseline_time * baseline_workers / w for w in workers]
    ax1.plot(workers, ideal_times, 'r--', linewidth=2, alpha=0.7,
             label='Ideal linear scaling')
    
    ax1.set_xlabel('Number of Workers', fontsize=12)
    ax1.set_ylabel('Total Execution Time (ms)', fontsize=12)
    ax1.set_title('Total Execution Time vs Number of Workers', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Add values on points
    for w, t in zip(workers, total_times):
        ax1.text(w, t, f'{t:.1f} ms', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    # Chart 2: Speedup and Efficiency
    speedups = [baseline_time / t for t in total_times]
    efficiencies = [s / (w / baseline_workers) * 100 for s, w in zip(speedups, workers)]
    
    line1 = ax2.plot(workers, speedups, marker='o', linewidth=3, markersize=10,
                     color='green', label='Speedup')
    
    line2 = ax2.plot(workers, efficiencies, marker='s', linewidth=3, markersize=10,
                     color='orange', label='Parallel Efficiency (%)')
    
    ax2.set_xlabel('Number of Workers', fontsize=12)
    ax2.set_ylabel('Speedup / Efficiency (%)', fontsize=12)
    ax2.set_title('Speedup and Parallel Efficiency', fontsize=14, fontweight='bold')
    
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax2.legend(lines, labels, loc='best', fontsize=10)
    
    ax2.grid(True, alpha=0.3)
    
    for w, s, e in zip(workers, speedups, efficiencies):
        ax2.text(w, s, f'{s:.2f}x', ha='center', va='bottom', fontsize=8, color='black', fontweight='bold')
        ax2.text(w, e, f'{e:.1f}%', ha='center', va='top', fontsize=8, color='black', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')

def generate_scaling_report(data_dict, all_tasks, task_time_percentages, output_file='scaling_report.txt'):
    workers = sorted(data_dict.keys())
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(f"Worker counts analyzed: {workers}\n")
        
        baseline_workers = workers[0]
        baseline_time = data_dict[baseline_workers]['time-ns'].sum() / 1e6
        final_workers = workers[-1]
        final_time = data_dict[final_workers]['time-ns'].sum() / 1e6
        
        f.write("\n" + "="*70 + "\n\n")
        f.write("DETAILED SCALING RESULTS\n")
        f.write("-" * 70 + "\n\n")
        
        f.write(f"{'Workers':<10} {'Time (ms)':<15} {'Speedup':<12} {'Efficiency (%)':<15}\n")
        f.write("-" * 70 + "\n")
        
        for w in workers:
            total_time = data_dict[w]['time-ns'].sum() / 1e6
            speedup = baseline_time / total_time
            eff = (speedup / (w / baseline_workers)) * 100
            f.write(f"{w:<10} {total_time:<15.2f} {speedup:<12.2f} {eff:<15.2f}\n")
        
        all_tasks_data = []
        for task_key, times in all_tasks.items():
            if baseline_workers not in times:
                continue
            
            task_id, task_name = task_key
            baseline_task_time = times[baseline_workers]
            time_percent = task_time_percentages.get(task_key, 0.0)
            
            speedups_by_workers = {}
            efficiencies_by_workers = []
            
            for w in workers:
                if w in times and times[w] > 0:
                    speedup = baseline_task_time / times[w]
                    speedups_by_workers[w] = speedup
                    
                    ideal_speedup_w = w / baseline_workers
                    efficiency_w = (speedup / ideal_speedup_w) * 100
                    efficiencies_by_workers.append(efficiency_w)
                else:
                    speedups_by_workers[w] = 0.0
            
            avg_efficiency = np.mean(efficiencies_by_workers) if efficiencies_by_workers else 0.0
            
            if final_workers in times and times[final_workers] > 0:
                final_speedup = baseline_task_time / times[final_workers]
            else:
                final_speedup = 0.0
            
            all_tasks_data.append((task_name, speedups_by_workers, avg_efficiency, final_speedup, task_id, time_percent))
        
        def write_task_table(tasks_subset, title, f):
            f.write(f"\n{'='*70}\n\n")
            f.write(f"{title}\n")
            f.write("-" * 70 + "\n\n")
            
            min_workers = workers[0]
            max_workers = workers[-1]
            time_pct_header = f"Time%({min_workers}w)"
            eff_header = f"Eff%({min_workers}w→{max_workers}w)"
            
            header = f"{'Task':<50} {time_pct_header:<12} {eff_header:<15}"
            for w in workers:
                header += f" {'Speedup ' + str(w) + 'w':<12}"
            f.write(header + "\n")
            
            separator_width = 50 + 12 + 15 + 12 * len(workers) + len(workers) + 2
            f.write("-" * separator_width + "\n")
            
            for task, speedups_dict, avg_eff, _, task_id, time_pct in tasks_subset:
                task_display = f"{task[:40]} (ID:{task_id})"
                line = f"{task_display[:50]:<50} {time_pct:<12.2f} {avg_eff:<15.2f}"
                for w in workers:
                    speedup_val = speedups_dict.get(w, 0.0)
                    line += f" {speedup_val:<12.2f}"
                f.write(line + "\n")
        
        # Top 20 best scaling tasks
        best_scaling = sorted(all_tasks_data, key=lambda x: x[3], reverse=True)
        write_task_table(best_scaling[:20], "TOP 20 BEST SCALING TASKS", f)
        
        # Top 20 worst scaling tasks
        # Filter tasks with >10ms at baseline
        significant_tasks = []
        for task_key, times in all_tasks.items():
            if baseline_workers in times and times[baseline_workers] > 10:
                task_id, task_name = task_key
                for task_data in all_tasks_data:
                    if task_data[4] == task_id:
                        significant_tasks.append(task_data)
                        break
        
        slow_scaling = sorted(significant_tasks, key=lambda x: x[3])
        write_task_table(slow_scaling[:20], 
                        f"TOP 20 WORST SCALING TASKS\n(Among tasks with >10ms at {baseline_workers} workers)", 
                        f)
        
        all_tasks_by_id = sorted(all_tasks_data, key=lambda x: x[4])
        write_task_table(all_tasks_by_id,
                        f"COMPLETE TASK LIST ({len(all_tasks_by_id)} tasks)\n(Sorted by task ID / order of appearance)",
                        f)
          
def main():
    parser = argparse.ArgumentParser(
        description='Analyze parallel scaling performance across different worker counts',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python analyze_scaling.py file-1.cali=1 file-2.cali=2 file-4.cali=4 file-8.cali=8
  python analyze_scaling.py hard-report-1.cali=1 hard-report-2.cali=2 -n 15
  python analyze_scaling.py run-1.cali=1 run-4.cali=4 run-8.cali=8 --output-prefix strong_scaling
  
Argument format:
  Each file must be specified as: filename=X
  where X is the number of workers
  Example: hard-report-4.cali=4
        """
    )
    
    parser.add_argument('files', nargs='+',
                        help='Caliper files with format: filename=X (where X is worker count)')
    parser.add_argument('-n', '--top-tasks', type=int, default=10,
                        help='Number of top tasks to display (default: 10)')
    parser.add_argument('--output-prefix', default='scaling',
                        help='Prefix for output files (default: scaling)')
    
    args = parser.parse_args()
    
    data_dict = load_all_files(args.files)
    
    if len(data_dict) < 2:
        print("\nError: At least 2 different worker counts are required for scaling analysis")
        sys.exit(1)
    
    print(f"\nLoaded {len(data_dict)} configurations")
    print(f"Worker counts: {sorted(data_dict.keys())}")
    
    all_tasks, task_time_percentages = get_task_scaling_data(data_dict)
    print(f"Analyzing {len(all_tasks)} unique tasks")
    
    plot_task_scaling(all_tasks, data_dict, top_n=args.top_tasks,
                        output_file=f'{args.output_prefix}_tasks.png')
    plot_total_scaling(data_dict, output_file=f'{args.output_prefix}_total.png')
    generate_scaling_report(data_dict, all_tasks, task_time_percentages,
                            output_file=f'{args.output_prefix}_report.txt')
    
    print("\nScaling analysis completed successfully!")
    print(f"\nGenerated files:")
    print(f"  - {args.output_prefix}_tasks.png (task scaling)")
    print(f"  - {args.output_prefix}_total.png (total time and efficiency)")
    print(f"  - {args.output_prefix}_report.txt (detailed report)")

if __name__ == "__main__":
    main()