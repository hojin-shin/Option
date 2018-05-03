import math
import random
import matplotlib.pyplot as plt
import time
t0 = time.time()
s, k = 100, 100
r, q, t, v = 0.02, 0.01, 0.1, 0.25
flag = "put"
d = 1 if flag.lower()=="call" else -1

ps = [] 
ss = 0
nsim = 100000
for i in range(nsim):
    e = random.normalvariate(0,1)
    st = s * math.exp((r-q-0.5*v*v)*t \
                      + v*math.sqrt(t)*e)
    ss += max(d*(st-k),0) * math.exp(-r*t)
    
    #ps.append(st)

#plt.hist(ps, bins=40)
price = ss / nsim
print("price = {0}".format(price))
t1 = time.time()
print("time = {0}".format(t1-t0))
