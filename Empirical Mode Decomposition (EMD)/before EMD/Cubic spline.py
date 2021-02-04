import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

x = np.arange(0.2, 2*np.pi+np.pi/4, 2*np.pi/8)    # x coordinates of points
y = np.sin(x)                                   # y coordinates of points
tck = interpolate.splrep(x, y, s=0)          # Given the set of data points (x[i], y[i]) Find the B-spline representation of 1-D curve.
xnew = np.arange(0.2, 2*np.pi, np.pi/50)  # creates x values to build cubic spline
ynew = interpolate.splev(xnew, tck, der=0)  #evaluates spline for x points on y coordinates given by tck


plt.figure()
plt.plot(x, y, 'x', xnew, ynew, xnew, np.sin(xnew), x, y, 'b')
plt.legend(['Linear', 'Cubic Spline', 'True']) #linear comes from x,y and x,y, b
                                               # cubic spline comes from xnew, ynew
                                               #tru comes from xnew, np.sin(xnew)
plt.axis([-1.05, 7.33, -1.05, 1.05])
plt.title('Cubic-spline interpolation')
plt.show()