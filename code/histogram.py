from math import *
import matplotlib.pyplot as plt
import sys

#map key:energy value:count
energycounter = {}
with open(sys.argv[1]) as file:
    for i, line in enumerate(file):
        e = int(line)
        if e in energycounter.keys():
            energycounter[e] = energycounter[e]+1
        else:
            energycounter[e] = 1

listofkeys = []
listofcounts = []
for key in energycounter.keys():
    listofkeys.append(int(key))
    listofcounts.append(energycounter[key])

#plt.plot(listofkeys,listofcounts)
plt.bar(listofkeys,listofcounts,1)
plt.xlabel('Energy')
plt.ylabel('Number of events')
plt.title('P(E) T=2.4 Matrix=20x20 MCcycles=100000')
#plt.xticks(listofkeys+1*0.5,'label')
plt.savefig('Histogram_10000_10.png')
plt.show()
