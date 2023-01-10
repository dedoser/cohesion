import numpy as np;
import matplotlib.pyplot as plt
import sklearn as skl

energy = []
a = []
with open('res/e.txt', 'r') as e:
    energy = e.read().split(' ')

with open('res/a.txt', 'r') as aFd:
    a = aFd.read().split(' ')

xpoints = np.array(a)
ypoints = np.array(energy)


plt.plot(xpoints, ypoints)
plt.gca().invert_yaxis()
plt.xlim = a[-1]

plt.show()
plt.close()
