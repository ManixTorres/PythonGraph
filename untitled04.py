from cProfile import label
import enum
from itertools import count
import numpy
from operator import sub
import pandas as pd
import matplotlib as plt
import seaborn as sns
plt.use('TkAgg')

plt.style.use('fivethirtyeight')
MBP = pd.read_csv('alwoffii.tab', sep='\t', engine = 'python') ### name of file
MBP.to_csv('alwoffii.csv', index=None)

##seq_length = 2231536

nucleotides = 3259224
nuclooooo = list(range(1, 3259225))
#MBP = MBP.sort_values(by=['blastN'])

temp = 0
#for i in seq_length:


subject_from = MBP['subject from']  ##length
subject_to =  MBP['subject t0']
huey = MBP['%match']
id = MBP['table id | sample_id | DNA length (query)']



reads = abs(subject_from - subject_to)
PercentMatch =  reads
pReads = (reads * 100)  

for j, d in enumerate(subject_from): 
    if int(subject_to[j]) < int(subject_from[j]):
        temp = subject_from[j]
        subject_from[j] = subject_to[j]
        subject_to[j] = temp

#lolFirst.update({"id": index, "from": data, "to": subject_to[index]})

counting = []
counting = [-1 for i in range(nucleotides)]
clr = [0 for i in range(nucleotides)]
loopCount = 0


for n in range(0, 10133): #10133  
    for i in range (subject_from[n], subject_to[n]):
        counting[i] += 1
        clr[i] = huey[n]
            
"""
for m, rr in enumerate(subject_from):
    if counting[m] == -1:
        counting[m] = 0
        
      
for q in range(0, 3259223):        
    if clr[q] == 0:
        clr[q] = 70   
"""

for p in range(0, 3259223):
    if counting[p] == 0 or -1:
        del counting[p]
        del clr[p]
        del nuclooooo[p]
        
"""             
#for n in range(0, 10133): 
        #clr[i] = huey[n]    
####################################################################################################### 



# Dot plot created using scatter plot

cm = plt.cm.get_cmap("RdYlGn")
scale_legend = plt.Normalize(min(clr) + 1, max(clr)) 
color_map = plt.cm.ScalarMappable(cmap="RdYlGn", norm=scale_legend)
color_map.set_array([])

ax = plt.scatter(nuclooooo, counting, cmap="jet", s=60) 
ax.figure.colorbar(color_map)
#ax.get_legend().remove()
#ax.grid(False) 



#plt.colorbar(clr)
#plt.figure(figsize=(15, 10))
#plt.ylabel("Overlapping Read Depth", size=12)
#plt.xlabel("Alwoffii", size=12)
#plt.gcf().text(0.83, 0.5, "BlastN Percent Identity", fontsize=14, rotation=90)
plt.title("", size=20)
plt.tight_layout()
plt.show()

"""




ax = plt.scatter(nuclooooo, counting, s=60, c=clr, cmap="RdYlGn", alpha=0.1)
cbar = plt.colorbar(orientation = "vertical")
cbar.set_label(label = "identity")
plt.ylabel("Overlapping Read Depth", size=12)
plt.xlabel("Alwoffii", size=12)


plt.show()



