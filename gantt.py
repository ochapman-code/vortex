import numpy as np
import matplotlib.pyplot as plt

file = "mpi_img/gantt_chart_data.txt"      # File name

with open(file, mode='r') as file:          # Read file
    fileContent = file.readlines()

n = len(fileContent[1].split())
            
data = np.empty([len(fileContent[1:]), n])
for i, line in enumerate(fileContent[1:]):
    for j, val in enumerate(line.split()):
        data[i,j] = float(val)              # Enter data into array

processors = range(np.shape(data)[0])
segments = n

data=np.array(data).T

y_pos = np.arange(len(processors))

fig = plt.figure(figsize=(5,4))     # Start figure
ax = fig.add_subplot(111)

colors = fileContent[0].replace("\n","")    # Get colours from file
patch_handles = []

# left alignment of data starts at zero
left = np.zeros(len(processors)) 
for i, d in enumerate(data):
    patch_handles.append(ax.barh(y_pos, d, 
      color=colors[i%len(colors)], align='center', 
      left=left))
    left += d

from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# Add legend
legend_elements = [Patch(facecolor='g', edgecolor='g', label='Computation'),
                    Patch(facecolor='r', edgecolor='r', label='Idle'),
                    Patch(facecolor='b', edgecolor='b', label='Message passing'),
                    Patch(facecolor='c', edgecolor='c', label='Saving')]

# Create the figure
ax.legend(handles=legend_elements, loc='upper right')

ax.set_yticks(y_pos)
ax.set_yticklabels(processors)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Processor')
plt.tight_layout()
plt.show()

