from math import *
import matplotlib.pyplot as plt
import sys

x = []
y = []
with open(sys.argv[1]) as file:
    for i, line in enumerate(file):
        if line.startswith('>'):
            continue
        line = line.rstrip()
        data = line.split()
        x.append(float(data[int(sys.argv[2])]))
        y.append(float(data[int(sys.argv[3])]))

plt.plot(x,y)
if len(sys.argv) > 4:
    plt.xlabel(sys.argv[4])
    plt.ylabel(sys.argv[5])
else:
    plt.xlabel('x')
    plt.ylabel('y')
plt.show()

