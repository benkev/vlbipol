import numpy as np
from scipy.stats import norm, chi2

alpha = np.array([0.01, 0.025, 0.05, 0.95, 0.975, 0.99])

for df in range(1, 7):
    print("%2d %4.1f %4.1f %4.1f %6.4f %7.5f %7.5f" %
          (df, *tuple(chi2.ppf(1-alpha, df))))
    # print("%2d %4.1f %4.1f %4.1f %6.4f %7.5f %7.5f" %
    #       (df, *tuple(chi2.isf(alpha, df))))
      
for df in range(6, 31):
    print("%2d %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f" %
          (df, *tuple(chi2.ppf(1-alpha, df))))
    # print("%2d %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f" %
    #       (df, *tuple(chi2.isf(alpha, df))))


              

