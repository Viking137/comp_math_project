import numpy as np
import matplotlib.pyplot as plt

X = np.loadtxt('X.txt')#, delimiter= '\t', dtype = np.float)
Y = np.loadtxt('Y.txt')#, delimiter= '\t', dtype = np.float)

plt.plot(np.log(X), np.log(Y), '.')
plt.grid()
plt.show()
