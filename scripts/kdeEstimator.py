import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import numpy as np
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV, LeaveOneOut
import pandas as pd
from scipy import stats
from scipy.stats import gaussian_kde
barTab=pd.read_csv('outs/barTable.csv')
cellTemp = np.log2(barTab.iloc[:,4].values.astype(int))
cellTemp[cellTemp == -np.inf] = 0
cellTemp = stats.zscore(cellTemp)
cellTemp_d=np.linspace(np.quantile(cellTemp, 0.001), np.quantile(cellTemp, 0.999), 1000)

gkde_fit = gaussian_kde(cellTemp)
gkdeEval = gkde_fit.evaluate(cellTemp_d)
sns.distplot(gkdeEval, rug=True, hist=False)

machInt= np.iinfo(np.int32).min
y=np.append(machInt, gkdeEval)
y_=np.diff(y) > 0 

plt.figure()
# plot histgram of sample
plt.hist(cellTemp, bins=20, normed=1)
# plot data generating density
plt.plot(cellTemp, stats.norm.pdf(cellTemp), color="r", label='DGP normal')
# plot estimated density
plt.plot(cellTemp_d, gkdeEval, label='kde', color="g")
plt.title('Kernel Density Estimation')
plt.legend()
#plt.show()