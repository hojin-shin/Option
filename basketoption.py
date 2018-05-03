import numpy as np
import datetime
from scipy.stats import mvn, norm
import dateutil.relativedelta as rd

def euroTwoAssetMinBasketCall(f1, f2, k, df, v1, v2, rho):
    stdDev1, stdDev2 = v1, v2
    variance = v1*v1 + v2*v2 - 2*rho*v1*v2
    stdDev = np.sqrt(variance)
    modRho1 = (rho * stdDev2 - stdDev1) / stdDev
    modRho2 = (rho * stdDev1 - stdDev2) / stdDev
    D1 = (np.log(f1/f2) + 0.5*variance) / stdDev
    if (k==0):
        gamma = 1.0
        alfa = norm.cdf(-D1)
        beta = norm.cdf(D1 - stdDev)
    else:
        D1_1 = (np.log(f1/k) + 0.5*v1*v1) / stdDev1
        D1_2 = (np.log(f2/k) + 0.5*v2*v2) / stdDev2
        lower = np.array([-10,-10])
        mean = np.array([0,0])
        cov0 = np.array([[1,rho],[rho,1]])
        cov1 = np.array([[1,modRho1],[modRho1,1]])
        cov2 = np.array([[1,modRho2],[modRho2,1]])
        alfa = mvn.mvnun(lower, np.array([D1_1, -D1]), mean, cov1)[0]
        beta = mvn.mvnun(lower, np.array([D1_2, D1 - stdDev]), mean, cov2)[0]
        gamma = mvn.mvnun(lower, np.array([D1_1 - stdDev1, D1_2 - stdDev2]), mean, cov0)[0]
    return df * (f1*alfa + f2*beta - k*gamma)    
    
def euroTwoAssetMinBasketPut(f1, f2, k, df, v1, v2, rho):    
    return k * df - euroTwoAssetMinBasketCall(f1, f2, 0.0, df, v1, v2, rho) + \
                    euroTwoAssetMinBasketCall(f1, f2, k, df, v1, v2, rho);

def worstOfPut(evaluationDate, maturityDate, k, s1, s2, r, q1, q2, sigma1, sigma2, rho):
    t = (maturityDate - evaluationDate).days / 365.0
    if t==0:
        return np.maximum((-1)*(np.minimum(s1,s2)-k),0)
    else:
        df = np.exp(-r*t)
        f1 = s1 / df * np.exp(-q1*t)
        f2 = s2 / df * np.exp(-q2*t)
        v1 = sigma1 * np.sqrt(t)
        v2 = sigma2 * np.sqrt(t)    
        return euroTwoAssetMinBasketPut(f1, f2, k, df, v1, v2, rho)

def bsbasket(f, evaluationDate, maturityDate, k, s1, s2, r, q1, q2, sigma1, sigma2, corr):    
    p = np.zeros((3,3))
    for i in range(-1,2):
        for j in range(-1,2):
            p[i+1][j+1] = \
             f(evaluationDate, maturityDate, k, s1+i, s2+j, r, q1, q2, sigma1, sigma2, corr)    
    npv = p[1][1]
    delta = np.array([p[2][1]-p[0][1], p[1][2]-p[1][0]]) / 2.0
    gamma = np.array([p[2][1]-2*npv+p[0][1], p[1][2]-2*npv+p[1][0]])
    xgamma = (p[2][2] + p[0][0] - p[2][0] - p[0][2]) / 4.0
    if (evaluationDate==maturityDate):
        theta = 0;
    else:
        theta = f(evaluationDate+rd.relativedelta(days=1), \
                  maturityDate, k, s1, s2, r, q1, q2, sigma1, sigma2, corr) - npv
    vega = np.array([f(evaluationDate, maturityDate, k, s1, s2, r, q1, q2, sigma1+0.01, sigma2, corr) \
                     - npv,
        f(evaluationDate, maturityDate, k, s1, s2, r, q1, q2, sigma1, sigma2+0.01, corr) - npv])    
    value = {"npv":npv, "delta":delta, "gamma":gamma, "xgamma":xgamma, "theta":theta, "vega":vega}
    return value
    
if __name__=="__main__":
    today = datetime.date(2018,5,3)
    maturity = datetime.date(2018,5,10)
    s1, s2, k = 100,100,100
    r = 0.02
    q1, q2 = 0.01, 0.01
    sigma1, sigma2 = 0.2, 0.2
    rho = 0.6
    value0 = bsbasket(worstOfPut, today, maturity, k, s1, s2, r, q1, q2, sigma1, sigma2, rho)
    print(value0)

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    
    ss1 = np.arange(90,111)
    ss2 = np.arange(90,111)
    X, Y = np.meshgrid(ss1, ss2)
    value = np.zeros((len(ss1), len(ss2)))
    for i,si in enumerate(ss1):
        for j,sj in enumerate(ss2):
            value[i][j] = (-1)*bsbasket(worstOfPut, today, maturity, k, si, sj, r, q1, q2, sigma1, sigma2, rho)['npv']
    fig = plt.figure(1)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, value, linewidth=1, rstride=1, cstride=1, cmap=cm.coolwarm)
    

    #cross-gamma effect
    '''
    ss1 = np.arange(90,111)
    ss2 = np.arange(90,111)
    X, Y = np.meshgrid(ss1, ss2)
    value = np.zeros((len(ss1), len(ss2)))
    dt = 5.0/252
    for i,si in enumerate(ss1):
        for j,sj in enumerate(ss2):
            value[i][j] = value0['xgamma']*((si-s1)*(sj-s2) - rho*sigma1*sigma2*s1*s2*dt)
    fig = plt.figure(1)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, value, linewidth=1, rstride=1, cstride=1, cmap=cm.coolwarm)   
    
   
    plt.figure(figsize=(5,5))
    CS = plt.contour(X, Y, value, 10)
    plt.clabel(CS)
    plt.title('Contour Plot')
    '''