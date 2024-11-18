from _imports import *

a = 2
b = 1

t = np.linspace(0, 2 * np.pi, 1000)
x = a * np.cos(t)
y = b * np.sin(t) 

plt.figure()
plt.plot(x, y)
plt.show()