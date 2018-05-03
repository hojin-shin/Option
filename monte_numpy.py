import math
import random
import matplotlib.pyplot as plt
import time
import numpy as np

t0 = time.time()
s, k = 100, 100
r, q, t, v = 0.02, 0.01, 0.1, 0.25
flag = "put"
d = 1 if flag.lower()=="call" else -1
nsim = 100000

e = np.random.randn(nsim)
st = s*np.exp((r-q-0.5*v*v)*t + v*np.sqrt(t)*e)
#plt.hist(st,bins=50)

payoff = np.maximum(d*(st-k),0) * np.exp(-r*t)
price = payoff.mean()

print("price = {0}".format(price))
print("time = {0}".format(time.time()-t0))
