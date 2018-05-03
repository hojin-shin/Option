
import numpy as np
import datetime as dt
import scipy.stats as sst

def bspv(today, matDate, k, s, r, q, sigma, optionType):
    t = (matDate - today).days / 365.0
    d1 = (np.log(s/k) + (r-q+0.5*sigma**2)*t) / (sigma*np.sqrt(t))
    d2 = d1 - sigma*np.sqrt(t)    
    cp = 1 if optionType.lower()=="call" else -1
    nd1 = sst.norm.cdf(cp*d1)
    nd2 = sst.norm.cdf(cp*d2)
    if today==matDate:
        npv = np.maximum(cp*(s-k),0)
        delta, gamma, theta, vega, rho = (np.zeros(s.shape),)*5
    else:
        npv = cp*(s*np.exp(-q*t)*nd1 - k*np.exp(-r*t)*nd2)  
        delta = cp * np.exp(-q*t)*nd1
        gamma = sst.norm.pdf(d1) * np.exp(-q*t) / (s * sigma * np.sqrt(t))
        theta = -s * sigma * sst.norm.pdf(d1) * np.exp(-q*t) \
                / (2 * np.sqrt(t)) + cp * q * s * nd1 * np.exp(-q * t) \
                - cp * r * k * np.exp(-r * t) * nd2
        vega = s * np.sqrt(t) * sst.norm.pdf(d1) * np.exp(-q*t)
        rho = cp * k * t * np.exp(-r * t) * nd2
   
    value = {"npv":npv, "delta":delta, "gamma":gamma, \
            "theta":theta, "vega":vega, "rho":rho}
            
    return value

if __name__=="__main__":
    evaluationDate = dt.date(2016,11,15)
    maturityDate = dt.date(2016,12,30)
    s = np.arange(90,111)
    k = 100
    r = 0.02
    q = 0.01
    sigma = 0.2
    value = bspv(evaluationDate,maturityDate,k,s,r,q,sigma,"call")
    
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2,3,figsize=(10,7))
    otype = [['npv','delta','gamma'],['rho','vega','theta']]
    for i in range(2):
        for j in range(3):
            ax[i][j].plot(s, value[otype[i][j]])
            ax[i][j].set_title(otype[i][j].upper())
plt.show()

