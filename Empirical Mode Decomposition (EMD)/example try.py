import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

# given values
xi = np.array([[0.2, 0.5, 0.7],
              [0.9,   1,   3]])
dat_list = []
yi = np.array([0.3, -0.1, 0.2, 0.1])
num = np.array([1, 2])

dat_list.append(xi)
dat_list.append(yi)

print(2*dat_list)

print(dat_list)


