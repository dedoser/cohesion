import numpy as np;
import matplotlib.pyplot as plt

energy = []
a = []
with open('res/e1.txt', 'r') as e:
    for word in e.read().split(' '):
        energy.append(float(word))

with open('res/a1.txt', 'r') as aFd:
    for word in aFd.read().split(' '):
        a.append(float(word))

a = a[20:]
energy = energy[20:]

xpoints = np.array(a)
ypoints = np.array(energy)

energy1 = []
a1 = []
with open('res/e2.txt', 'r') as e:
    for word in e.read().split(' '):
        energy1.append(float(word))

with open('res/a2.txt', 'r') as aFd:
    for word in aFd.read().split(' '):
        a1.append(float(word))
energy1 = energy1[20:]
a1 = a1[20:]

xpoints1 = np.array(a1)
ypoints1 = np.array(energy1)

energy2 = []
a2 = []
with open('res/e3.txt', 'r') as e:
    for word in e.read().split(' '):
        energy2.append(float(word))

with open('res/a3.txt', 'r') as aFd:
    for word in aFd.read().split(' '):
        a2.append(float(word))

a2 = a2[15:]
energy2 = energy2[15:]

xpoints2 = np.array(a2)
ypoints2 = np.array(energy2)

plt.plot(xpoints, ypoints, label = "B-B")
plt.plot(xpoints1, ypoints1, label = "A-B")
plt.plot(xpoints2, ypoints2, label = "A-A")

# plt.gca().invert_yaxis()

plt.legend()
plt.show()
plt.close()
