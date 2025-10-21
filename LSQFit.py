import numpy as np
from numpy.linalg import inv
from matplotlib import pyplot as plt

xmin=1.0
xmax=20.0
npoints=12
sigma=0.2
lx=np.zeros(npoints)
ly=np.zeros(npoints)
ley=np.zeros(npoints)
pars=[0.5,1.3,0.5]

from math import log
def f(x,par):
    return par[0]+par[1]*log(x)+par[2]*log(x)*log(x)

from random import gauss
def getX(x):  # x = array-like
    step=(xmax-xmin)/npoints
    for i in range(npoints):
        x[i]=xmin+i*step

def getY(x,y,ey):  # x,y,ey = array-like
    for i in range(npoints):
        y[i]=f(x[i],pars)+gauss(0,sigma)
        ey[i]=sigma

# get a random sampling of the (x,y) data points, rerun to generate different data sets for the plot below
getX(lx)
getY(lx,ly,ley)

fig, ax = plt.subplots()
ax.errorbar(lx, ly, yerr=ley)
ax.set_title("Pseudoexperiment")
fig.show()


# *** modify and add your code here ***
#chi2 and for nr loop designed with chatgpt
def bestfit(x,y,err):
  nPnts = len(lx)
  nPar  = 3

  A=np.matrix(np.zeros((nPnts, nPar)))
  #Fill the A matrix
  for nr in range(nPnts):
    A[nr,0] = 1
    A[nr,1] = lx[nr]
    A[nr,2] = lx[nr]*lx[nr]

  for i in range(nPnts):
    A[i] = A[i] / ley[i]
  yw = (ly/ley).reshape(nPnts,1)
  theta =  inv(np.transpose(A).dot(A)).dot(np.transpose(A)).dot(yw)

  a = theta[0,0]
  b = theta[1,0]
  c = theta[2,0]

  y_fit = a + b*np.log(x)+c*np.log(x)*np.log(x)
  chi2 = np.sum(((y - y_fit)/err)**2)
  chi2_reduced = chi2 / (npoints - 3)  
  return a,b,c,chi2_reduced
# perform many least squares fits on different pseudo experiments here
# fill histograms w/ required data
nexperiments = 1000

par_a = np.zeros(nexperiments)   # simple placeholders for making the plot example
par_b = np.zeros(nexperiments)   # these need to be filled using results from your fits
par_c = np.zeros(nexperiments)
chi2_reduced = np.zeros(nexperiments)


for i in range(nexperiments):
    getX(lx)
    getY(lx, ly, ley)
    par_a[i],par_b[i],par_c[i],chi2_reduced[i] = bestfit(lx,ly,ley)
    
    
fig, axs = plt.subplots(2, 2)
plt.tight_layout()

# careful, the automated binning may not be optimal for displaying your results!
axs[0, 0].hist2d(par_a, par_b)
axs[0, 0].set_title('Parameter b vs a')

axs[0, 1].hist2d(par_a, par_c)
axs[0, 1].set_title('Parameter c vs a')

axs[1, 0].hist2d(par_b, par_c)
axs[1, 0].set_title('Parameter c vs b')

axs[1, 1].hist(chi2_reduced)
axs[1, 1].set_title('Reduce chi^2 distribution')

fig.show()

# **************************************


input("hit Enter to exit")
