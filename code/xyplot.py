from math import *
import matplotlib.pyplot as plt
import sys

temp = []
energy = []
with open(sys.argv[1]) as file:
    for i, line in enumerate(file):
        if line.startswith('>'):
            continue
        line = line.rstrip()
        data = line.split()
        temp.append(float(data[0]))
        energy.append(float(data[1]))

plt.plot(temp,energy)
plt.xlabel('Temperature')
plt.ylabel('Energy')
plt.axis([1.6,2.6,-2.5,0])
plt.show()
