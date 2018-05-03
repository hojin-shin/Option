from blackscholes import bspv 
import datetime

s0, k, r, q, vol = 100, 100, 0.02, 0.01, 0.2
today = datetime.date(2016,11,15)
maturity = datetime.date(2016,12,15)
optiontype = "put"

value0 = bspv(today, maturity, k, s0, r, q, vol, optiontype)
acc = -value0['npv']
stock = -value0['delta']*s0
bond = acc - stock

print(acc, stock, bond)


import numpy as np
numSim = 1000
days = (maturity - today).days
dt = 1/365.0
e = np.random.randn(numSim, days)
s = np.zeros((numSim, days+1))
s[:,0] = s0

realizedVol = 0.2
drift = (r-q-0.5*realizedVol**2)*dt
for i in range(days) :
    s[:,i+1] = s[:,i] * np.exp(drift + realizedVol*np.sqrt(dt)*e[:,i])

'''    
import matplotlib.pyplot as plt
plt.plot(s.T)
'''

import dateutil.relativedelta as rd
d = today
c = 0
delta = np.ones(numSim) * (-value0['delta'])
final = np.zeros(numSim) 

while True:
    c += 1
    d = d + rd.relativedelta(days = 1)
    print(d)
    acc += (bond*(np.exp(r*dt)-1) + delta*(s[:,c]-s[:,c-1]))
    if d == maturity:
        for i in range(numSim):
            value = bspv(d, maturity, k, s[i,c], r, q, vol, optiontype)
            final[i] = -value['npv']
    break 

    for i in range(numSim):
        value = bspv(d, maturity, k, s[i,c], r, q, vol, optiontype)
        delta[i] = -value['delta']
    bond = acc - delta*s[:,c]

import matplotlib.pyplot as plt
plt.plot(s[:,-1],acc,'.')

he = acc-final
plt.plot(s[:,-1],he,'.')
plt.show()


