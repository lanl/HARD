"""
@Shihab
Analyzes results from strong/weak scaling runs.
"""

import os
import matplotlib.pyplot as plt 
import numpy as np 


#-----Multiple node, each 36 ranks, 512^3 box per node-------
MESH = 512
nodes = [1, 2, 4, 8, 16]
N = len(nodes)
mesh = [[9, 9, 9], [10, 9, 9], [10, 10, 9], [10, 10, 10],
        [11, 10, 10]]
secs = [[79, 109, 82, 86, 90], [78, 80, 82, 115, 118]]

secs = [(secs[0][i]+secs[1][i])/2 for i in range(N)]
print(secs)

fig, ax1 = plt.subplots()
ax1.set_title(f"{MESH}x{MESH}x{MESH} Mesh, Each Node has 36 ranks")
ax1.set_xlabel("# Nodes")
filename = 'weak.png'

rank = nodes
#----------------------------------------
ideal = [secs[0] for r in range(0,N)]
speedup = [1] + [secs[0]/s for s in secs[1:]]
efficiency = [1]+[secs[0]/s for s in secs[1:]]

#---------------------------------------------


ax1.set_ylabel("Seconds")
ax1.set_xscale("log")
ax1.set_ylim(np.min(secs)-10, np.max(secs)+10)
#ax1.set_yscale("log")

# Plot runtime and ideal on the primary y-axis
ax1.plot(rank, secs, label='Runtime', color='blue')
ax1.scatter(rank, secs, color='blue')  # Scatter plot of the same data
ax1.plot(rank, ideal, '--', label='Ideal Runtime', color='green')

# Set custom ticks
ax1.set_xticks(np.unique(rank))
ax1.set_xticklabels([str(r) for r in np.unique(rank)])

primary_y_ticks = np.geomspace(np.min(secs), np.max(secs), num=5) 
ax1.set_yticks(primary_y_ticks)
ax1.set_yticklabels([f"{y:.2f}" for y in primary_y_ticks])

# Create a secondary y-axis for efficiency
ax2 = ax1.twinx()
ax2.set_ylabel('Efficiency', color='y')  # Label the second y-axis
ax2.plot(rank, efficiency, 'o-', label='Efficiency', color='y')
ax2.tick_params(axis='y', labelcolor='y')

# Adding legend that combines both axes
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc='lower left')

# Save the figure
plt.savefig(filename)
plt.show()