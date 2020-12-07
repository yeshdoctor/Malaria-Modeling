import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
X0 = 0 #initial infected RBC count
I0 = 1e-2 #initial immune cell conc
A0 = 0 #initial antibody conc
M0 = 3e13 
T0 = 0
F0 = 0
tspan = np.linspace(0,20,100001)
k0 = 9e-5
k1 = 2e-9
k2 = 6e-4
ax = 0.025
k3 = 10e-8
r = 16
k4 = 8.5e-4
ay = 48
k5 = 10e-8
lamI = 10
lamx = 0.05
lamy = 0.05
ai = 0.05
k6 = 2000
k7 = 1500
nu = 0.6
aa = 5e-10
k8 = 1500 
c = -(2.5)
s2 = 2
f = 5.2e-15

tspan2 = np.linspace(0,48,10000)
for i in range (4):
    tspan2 = np.append(tspan2,tspan2)
tspan2 = tspan2[:85000]
tspan2 = np.append(np.zeros(15000),tspan2)

#%% Liver
spori = 16
tspanliv = np.linspace(0,2,10000)
mern = spori*np.exp(4.95174377627*tspanliv)
Y0 = np.max(mern) # initial merozoite #
y0 = [X0,Y0,I0,A0,M0,F0]

#%% main function
def s(t):
  return 0.5*(np.sin(np.pi*t/2)**2)
def s2(t):
    return(t%2)
def s4(t):
    if(s2(t)*24<24.058):
        s3=.0046*s2(t)*24+.0215
    else:
        s3=.1422*np.log(s2(t)*24-24)+.5382
    return s3
def iRBCcount(y0,t):
  ########################## Establishing local var names
  x = y0[0] 
  y = y0[1]
  i = y0[2]
  a = y0[3]
  m = y0[4]
  f = y0[5]
  ########################### all kinetic rate constants, meanings are defined in the paper. Taken directly from the paper. 
  dXdt = (k1/(1+k2*a))*y*m - s(t)*x-ax*x-k3*i*x-f
  dYdt = (r/(1+k4*i))*s(t)*x-ay*y-k5*i*y - (k1/(1+k2*a))*y*m
  dIdt = lamI+(((lamx*x)/(k6 + x)) + ((lamy*y)/(k7 + y)))*i - ai*i
  dAdt = nu*i*(y/(k8 + y)) - aa*a*y
  dMdt = (k0*x-.0083*m-k1/(1+k2*a)*y*m)*50#assuming 50 days age to avoid PDEs listed in paper
  dFdt = (0.16*s2(t)**.72)*(dXdt)*.047
  return [dXdt,dYdt,dIdt,dAdt,dMdt,dFdt]
vals = odeint(iRBCcount, y0, tspan)
lamI = 100
X0 = 1e6 #https://microbiologynotes.com/differences-between-primary-and-secondary-immune-response/
vals2 = odeint(iRBCcount, y0, tspan)
#%% sequestration
seq = np.zeros(100001)
circ = np.zeros(100001)
s3 = np.zeros(100001)
x = vals[:,0]
    
for i in range (0,100000):
    if(tspan2[i]<24.058):
        s3[i]=.0046*tspan2[i]+.0215
    else:
        s3[i]=.1422*np.log(tspan2[i]-24)+.5382
    seq[i] = s3[i]*x[i]
    circ[i] = (1-s3[i])*x[i]

# %% HRP2
def hrp2(x):
    hrp2 = np.zeros(len(tspan))
    nT = 0
    nmax = 0
    for i in range(len(tspan)):
        if tspan[i]%2 == 0 and nT != 0:#assuming cycle period of 2 days
            nT+=1
            hrp2[i] += hrp2[i-1]+f*x[i]
            nmax = i
        else:
            hrp2[i] = hrp2[nmax]*np.exp(-.18873*(tspan[i-nmax]))
        if nT == 0:
            nT = 1
    
    return hrp2/5 #per 5L blood
#%% Simulation
"""
nums = np.linspace(6e-4,1,20) #k2 range
nums2 = np.linspace(0,100,20) #lamI range
X,Y=np.meshgrid(nums/6e-4,nums2/10)

rows, cols = (20, 20) 
Z1 = [[0 for i in range(cols)] for j in range(rows)]
Z2 = [[0 for i in range(cols)] for j in range(rows)]
Z3 = [[0 for i in range(cols)] for j in range(rows)]

for i in range(0,20):
    k2 = nums[i]
    for j in range(0,20):
        lamI = nums2[j]
        vals3 = odeint(iRBCcount, y0, tspan)
        #Z1[i][j] = 100-((3e13-max(vals3[:,0]))/3e13)*100
        #Z2[i][j] = tspan[np.argmax(vals3[:,0])]
   \
        occurences = np.where(hrp2(vals3[:,0])*1000000>= 100)
        if(len(occurences[0]) == 0):
            Z3[i][j] = 0
        else:
            d1 = occurences[0][len(occurences[0])-1]
            d2 = occurences[0][0]
            Z3[i][j] = tspan[d1]-tspan[d2]


fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.plot_surface(X,Y,np.asarray(Z1),cmap='inferno')
#ax = fig.add_subplot(111, projection='3d')
#ax.plot_surface(X,Y,np.asarray(Z2),cmap='inferno')
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X,Y,np.asarray(Z3),cmap='inferno')
for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(18)
for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(18)
for t in ax.zaxis.get_major_ticks(): t.label.set_fontsize(18)
"""

#%% plotting


"""
fig,axs = plt.subplots(3,1, figsize=(20,22),sharex=True)
fig.subplots_adjust(hspace=0.075)
sns.set()
sns.set_style("whitegrid")
sns.set_context("talk")
axes_style = {'linewidth':5}
axs[0].plot(tspan, vals[:,1], sns.xkcd_rgb["pale red"],label = "first infection", linewidth = 5)
axs[0].plot(tspan, vals2[:,1],  label = "second infection", linewidth = 5)
axs[1].plot(tspan, vals[:,2], sns.xkcd_rgb["pale red"], label = "first infection", linewidth = 5)
axs[1].plot(tspan, vals2[:,2],  label = "second infection", linewidth = 5)
axs[2].plot(tspan, vals[:,3], sns.xkcd_rgb["pale red"], label = "first infection", linewidth = 5)
axs[2].plot(tspan, vals2[:,3],  label = "second infection", linewidth = 5)
axs[2].tick_params(axis='x', which='major', labelsize=32)
axs[2].tick_params(axis='y', which='major', labelsize=24)
axs[1].tick_params(axis='y', which='major', labelsize=24)
axs[0].tick_params(axis='y', which='major', labelsize=24)
plt.savefig("F2PBL2.png",dpi=500)
"""


"""
sns.set()
sns.set_style("whitegrid")
sns.set_context("talk")
plt.figure(figsize=(16,8))
plt.plot(tspan, hrp2(vals[:,0])*1000000, sns.xkcd_rgb["pale red"], label = "first infection", linewidth = 5)#1000000 ng/mL per g/L
plt.plot(tspan, hrp2(vals2[:,0])*1000000,  label = "second infection", linewidth = 5)
plt.plot(tspan, np.full(len(tspan),100),'k--', label = "Detection Threshold",linewidth = 5)
plt.xticks(fontsize=32)
plt.yticks(fontsize=24)
plt.savefig("F4PBL2.png",dpi=500)
"""


"""
sns.set()
sns.set_style("white")
sns.set_context("talk")
f, (ax, ax2) = plt.subplots(2, 1, sharex=True,figsize=(16,8))
f.subplots_adjust(hspace=0.3)
ax.plot(tspan, vals[:,0]/1e11, sns.xkcd_rgb["deep red"], label = "Infected - first infection", linewidth = 5)
ax.plot(tspan, vals2[:,0]/1e11, sns.xkcd_rgb["deep blue"], label = "Infected - second infection", linewidth = 5)
ax.plot(tspan, (3e13 - vals[:,0])/1e13, sns.xkcd_rgb["bright red"], label = "Healthy - first infection", linewidth = 5)
ax.plot(tspan, (3e13 - vals2[:,0])/1e13,sns.xkcd_rgb["bright sky blue"],   label = "Healthy - second infection", linewidth = 5)
ax.plot(tspan, np.full(np.array(100001),3.0),sns.xkcd_rgb["aqua green"],   label = "Healthy Person", linewidth = 5)

ax2.plot(tspan, vals[:,0]/1e11, sns.xkcd_rgb["bright red"],linestyle = "solid", label = "Infected - first infection", linewidth = 5)
ax2.plot(tspan, vals2[:,0]/1e11, sns.xkcd_rgb["bright sky blue"], linestyle = "solid",label = "Infected - second infection", linewidth = 5)
ax2.plot(tspan, (3e13 - vals[:,0])/1e13, sns.xkcd_rgb["bright red"], linestyle = "solid", label = "Healthy - first infection", linewidth = 5)
ax2.plot(tspan, (3e13 - vals2[:,0])/1e13,sns.xkcd_rgb["bright sky blue"], linestyle = "solid",   label = "Healthy - second infection", linewidth = 5)
ax2.plot(tspan, np.zeros(100001),sns.xkcd_rgb["aqua green"], linestyle = "solid",  label = "Healthy Person", linewidth = 5)

ax.set_ylim(2.985, 3.001)  # outliers only
ax2.set_ylim(-0.1, 1.6)  # most of the data

ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.spines['left'].set_visible(True)
ax2.spines['left'].set_visible(True)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()

d = .0075  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

ax2.tick_params(axis='x', which='major', labelsize=26)
ax2.tick_params(axis='y', which='major', labelsize=20)
ax.tick_params(axis='y', which='major', labelsize=20)
plt.savefig("F1PBL2.png",dpi=300)
"""



fig,axs = plt.subplots(2,1)
fig.subplots_adjust(hspace=0.5)
sns.set()
sns.set_style("white")
sns.set_context("talk")
axes_style = {'linewidth':5}

axs[0].plot(tspan, seq, sns.xkcd_rgb["pale red"], label = "Sequestered", linewidth = 5,color=sns.xkcd_rgb["coral"])
axs[0].plot(tspan, circ,  label = "Circulating", linewidth = 5,color=sns.xkcd_rgb["blue green"])
axs[0].plot(tspan, circ+seq,  label = "Total iRBCs",color='k', linewidth = 5)
axs[1].plot(tspan, vals[:,5]*(1-s3), sns.xkcd_rgb["bright red"], label = "first infection", linewidth = 5)
axs[1].plot(tspan, vals2[:,5]*(1-s3),  label = "second infection", linewidth = 5,color=sns.xkcd_rgb["bright sky blue"] )



axs[1].tick_params(axis='x', which='major', labelsize=34)
axs[0].tick_params(axis='x', which='major', labelsize=34)
axs[1].tick_params(axis='y', which='major', labelsize=26)
axs[0].tick_params(axis='y', which='major', labelsize=26)
#plt.savefig("F3PBL2.png",dpi=500)



