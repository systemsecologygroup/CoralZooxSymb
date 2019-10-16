# -*- coding: utf-8 -*-
from __future__ import division

from scipy import exp, linspace, array, zeros, e, sqrt, mean, ones
from numpy import amin, amax, meshgrid
from scipy.integrate import odeint
from matplotlib import pyplot as plt
import matplotlib

matplotlib.rc('font',**{'family':'Times', 'sans-serif':['Times']})
matplotlib.rc('text', usetex = True)

# Original model equation
"""
fig = plt.figure()
plt.text(0.05, 0.15, r"$\displaystyle\frac{dH}{dt} = F(U,\kappa\,E, H)\,H$"
                      +"\n"+ r"$S = \Gamma\,H, \quad  0 \leq \,\Gamma \leq \,1$"
                      +"\n"+r"$\displaystyle\frac{dU}{dt} = N \displaystyle\frac{\partial F}{\partial U}$"
                      +"\n"+r"$F(U, \kappa\,E, H) = G_H\,\kappa\,E\,\left(1-\displaystyle\frac{H}{K_H}\right) - \mu - M_H $"
                      +"\n"+r"$\kappa\,E,$ where $\kappa(\Gamma) = \displaystyle\frac{\Gamma}{\Gamma_h + \Gamma}$ and $E(U) = \left(1-e^{-\beta U}\right)$"
                      +"\n"+r"$\mu(U, \Gamma) = \alpha\,U\,\left(\displaystyle\frac{\Gamma}{1- \Gamma}\right)$", fontsize=20)
fig.savefig('Desktop/Equations.pdf')#, transparent=True) # Change this line appropriately when saving file in another directory
plt.show()
"""
# Analysis of the NON-dimensional model 
# Parameters
A = 0.1 # parameter of cost of symbiosis
M = 0.15 # corals's natural mortality
l = 1 # this is not important here but it is 1 for consistency with the model 
SIGMA = 1 # this is the speed of the adaptive process in the model

delta = 100000  # avoid negative investements

fsize= 24 # fontsize for figures

def Ydensity(gamma1, X):
    return gamma1*X
    
gammaHalfOpt = 0.3  # this should be gamma_h half saturation 

m = 1
def FitnessX(U, X, gamma1):
    #Y = Ydensity(gamma1, X)
    E = 1 - exp(-U) # energy received from symbiosis
    kappa = gamma1**l/(gammaHalfOpt**l + gamma1**l)
    print kappa
    Benefit = kappa*E*(1-X)
    cost_gam = gamma1**m/(1-gamma1**m)
    Cost =  (1-exp(-delta*U))*A*U*(cost_gam) + M
    Fitness = (Benefit - Cost)
    return Fitness
    
# plotting fitness function as a function of U for a fixed value of X = coral biomass (which is H in ms1)
# Ms1 Figure 3
"""
U = linspace(0, 5, 100)
X = 0.2
gammaList = array([0.05, 0.25, 0.40, 0.7])
#color_list=array([(0.8, 0.8, 0.8),(0.55, 0.55, 0.55), (0.35, 0.35, 0.35), (0, 0, 0)]) # grey scale
color_list = ((0.482, 0.867, 0.867), (0.259, 0.608, 0.608), (0.145, 0.471, 0.471), (0.031, 0.376, 0.376))  # blue scale

sub = plt.subplot(2, 3, 1)
#plt.title("Host fitness at $H_{*}$=%.2f \n"%X, fontsize=30)
for i in xrange(len(gammaList)):
    gamma1 = gammaList[i]
    FIT = FitnessX(U, X, gamma1)
    sub.plot(U, FIT, label="$\gamma$ = %.2f"%gamma1, linewidth=4, color=color_list[i])
        
plt.plot(U, 0*U,"--", color= "k", linewidth = 3)
plt.legend(loc="best", fontsize=20)
plt.xticks((0, 1, 2, 3, 4, 5), ("$0$", "$1$", "$2$", "$3$", "$4$", "$5$"),fontsize=23)
plt.yticks((-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2), ("$-1.0$", "$-0.8$", "$-0.6$","$-0.4$",  "$-0.2$", "$0.0$", "$0.2$"), fontsize=23)
plt.xlabel("Trait value $u$", fontsize = fsize)
plt.ylabel("Holobiont fitness $f$", fontsize = fsize)
plt.xlim((0, max(U)))
plt.show()
"""
# plotting fitness in 2 dimension (MS1 Figure 4 & 5, need to uncomment relevent pieces of codes)

gammaList = array([0.08, 0.35, 0.65])
#title_list = array(["(a)", "(b)", "(c)"])
title_list = array(["A", "B", "C"])

UList, XList = meshgrid(linspace(0, 5, 300), linspace(0, 1, 200))

UList2, XList2 = meshgrid(linspace(0, 5, 20), linspace(0, 1, 20)) # for the vector field
NullLevels=array([-0.0000005, 0.0000005])


#sub = plt.subplot(2, 3, 4)
# white box behind the legends
#sub.text(0.01, 1.125, ", $\quad$ $\quad$ $\quad$ $\quad$ $\quad$$\quad$$\quad$ $\quad$  $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$$\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$$\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ $\quad$ ju s t \n, \n , \n ,"
#            , color="white", bbox=dict(facecolor = "white", alpha = 1), fontsize = 17)
            
"""
for i in xrange(len(gammaList)):
    gamma1 = gammaList[i]
    FITNESS = FitnessX(UList, XList, gamma1)
    
    levels = linspace(-1, 1, 500)
    sub= plt.subplot(2,3, i+1+3)
    if i+1 == 1:
        plt.text(-1, 0.9, title_list[i], fontsize = 25)
    else:
        plt.text(-0.6, 0.9, title_list[i], fontsize = 25)
        
    kappa = gamma1**l/(gammaHalfOpt**l + gamma1**l) 
    cost_gam = gamma1**m/(1-gamma1**m)
    FitGrad = -delta*A*UList*(cost_gam)*exp(-delta*UList) - A*(cost_gam)*(1 - exp(-delta*UList)) + kappa*(-XList + 1)*(exp(-UList))
    
    
    line1= plt.contour(UList, XList, FITNESS, NullLevels, colors=((0.145, 0.471, 0.471),), linewidths=5)
    line1.collections[1].set_label("Host nullcline $(dh = \,0)$ ")
    line2= plt.contour(UList,XList, FitGrad , NullLevels,colors=((0.482, 0.867, 0.867), ), linewidths=5)
    line2.collections[1].set_label("Trait nullcline $(du = \,0)$ ")
    
    artistLine1, labelsLine1_def = line1.legend_elements()
    artistLine2, labelsLine2_def = line2.legend_elements()
    labelsLine1 = [r"Host nullcline $\big{(dh = \,0)}$"]
    labelsLine2 = ["Trait nullcline $(du = \,0)$"]
    
    plt.xlabel("Trait value $u$", fontsize=fsize)
    if i+1==1:
        plt.ylabel("Host $h$", fontsize = fsize)
        plt.legend(borderpad=1.3, labelspacing =0.5, frameon = False, handlelength =3, handleheight=1.5, loc=(-0.05, 1.05), fontsize= 19)
        
    # Vector Field
   
    FITNESS2 = FitnessX(UList2, XList2, gamma1)
    cost_gam = gamma1**m/(1-gamma1**m)
    FitGrad2 = -delta*A*UList2*(cost_gam)*exp(-delta*UList2) - A*(cost_gam)*(1 - exp(-delta*UList2)) + kappa*(-XList2 + 1)*(exp(-UList2))
    #plt.quiver(UList2, XList2, FitGrad2/sqrt(FitGrad2**2 + FITNESS2**2),  FITNESS2/sqrt(FitGrad2**2 + FITNESS2**2), color="grey")#color=(0.031, 0.376, 0.376)) # vector field
    
   
    plt.xticks((0, 1, 2, 3, 4, 5), ("$0$", "$1$", "$2$", "$3$", "$4$", "$5$"), fontsize=fsize)
    if i+1==1:
        plt.yticks((0.0, 0.2, 0.4, 0.6, 0.8, 1.0), ("$0.0$", "$0.2$","$0.4$", "$0.6$", "$0.8$", "$1.0$"), fontsize=fsize) 
    else:
        plt.yticks((0.0, 0.2, 0.4, 0.6, 0.8, 1.0), (" ", " "," ", " ", " ", " "), fontsize=fsize)     

    # Comment from this part to the legend and uncoment plt.quiver (line 133 above) to get the vector field
 
    cs = plt.contourf(UList, XList, FITNESS, array([-20, 0, 20]),colors=((0.031, 0.376, 0.376), (0.259, 0.608, 0.608)),#colors=((0.75, 0.75, 0.75), (0.9, 0.9, 0.9)),
                  hatches=[' ', " "],
                  extend='lower'
                  )# File can only be saved as eps format
    # create a legend for the contour set
    if i == 1: # only put legend once
        artists1, labels_X= cs.legend_elements()
        labels1 = [r"Host Fitness $<0$ $(dh < \,0)$", "Host Fitness $>0$ $(dh < \,0)$"] # labels_X gives default labels that's not what is needed 
        plt.legend(artists1, labels1, borderpad=1.3, handlelength =3, handleheight=1.5, frameon = False, loc=(-0.05, 1.05), fontsize=19)
    cs2 = plt.contourf(UList, XList, FitGrad, array([-20, 0, 20]), colors='none',
                  hatches=['/', '\\'],
                  extend='lower'
                  )# File can only be saved as eps format
                  
    # create a legend for the contour set
    if i==2: # only put legend once
        artists2, labels_X= cs2.legend_elements()
        labels2 = [r"Fitness gradient $<0$ $(du < \,0)$", "Fitness gradient $>0$ $(du > \,0)$"] # labels_X gives default labels that's not what is needed 
        plt.legend(artists2, labels2, borderpad=1.3, handlelength =3, handleheight=1.5, frameon = False, loc=(-0.05, 1.05), fontsize=19)
  

#plt.legend([artistLine1[1]] + [artistLine2[1]] + artists1 + artists2, labelsLine1 + labelsLine2 + labels1 +labels2,columnspacing=16.75, handletextpad= 2, handlelength = 3, frameon = True,  handleheight =1.5,loc= (-2.4, 1.05), ncol=3, fontsize=15)
fig.savefig('Desktop/2DFitn.pdf', transparent=True) # Change this line appropriately when saving file in another directory

plt.show()
"""

# System dynamics
m=1
def Gradient(X, U, A, gamma1):
    kappa = gamma1**l/(gammaHalfOpt**l + gamma1**l) 
    cost_gam = gamma1**m/(1-gamma1**m)
    Grad = -delta*A*U*(cost_gam)*exp(-delta*U) - A*(cost_gam)*(1 - exp(-delta*U)) + kappa*(-X + 1)*(exp(-U))
    return Grad
    
def System(Initial, T, A, gamma1):
    dSystem = zeros(len(Initial))
    X = Initial[0]
    U = Initial[1]
    E = 1 - exp(-U)
    kappa = gamma1**l/(gammaHalfOpt**l + gamma1**l) # Figure
    Benefit = kappa*E*(1-X)
    cost_gam = gamma1**m/(1-gamma1**m)
    Cost =  (1-exp(-delta*U))*A*U*(cost_gam) + M
    Fitness = (Benefit - Cost)
    dSystem[0] = Fitness*X
    dSystem[1] = SIGMA*Gradient(X, U, A, gamma1)
    return dSystem
     
Initial = array([0.8, 2]) 

LineType = array(["--", "-", "--", "--","--", "--"])

Timemax = 1000
TimeList = linspace(0., Timemax, 10000)



# Temporal variation (Ms1 Figure 6)
"""
#color_list=array([(0.8, 0.8, 0.8),(0.55, 0.55, 0.55), (0.25, 0.25, 0.25)])
color_list = ((0.482, 0.867, 0.867), (0.259, 0.608, 0.608), (0.145, 0.471, 0.471), (0.031, 0.376, 0.376))  # blue scale
color_list = ((0.482+0.5, 0.867, 0.867), (0.259+0.5, 0.608-0.1, 0.608-0.1), (0.145+0.5, 0.471-0.1, 0.471), (0.031+0.5, 0.376-0.2, 0.376-0.2))  # red scale
color_list = ((0.898, 0.60, 0.678), (0.857, 0.276, 0.278), (0.333, 0, 0)) # for low, mid and high gamma

label_list=array(["Low", "Mid", "High"])

sub1 = plt.subplot(2, 3, 1)
plt.text(-150, 4.5, title_list[0], fontsize = 25)
plt.xlabel(r"Model time $\tau$", fontsize=fsize)
plt.ylabel("Trait value $u$", fontsize=fsize)
plt.xticks((0, 200, 400, 600, 800, 1000), ("$0$", "$200$", "$400$", "$600$", "$800$", "$1000$"), fontsize=22)
plt.yticks((0, 1, 2, 3, 4, 5), ("$0$", "$1$", "$2$", "$3$", "$4$", "$5$"), fontsize=22)
plt.ylim((0, 5))
plt.xlim((0, Timemax))

sub2 = plt.subplot(2, 3, 2)
plt.text(-165, 0.9, title_list[1], fontsize = 25)
plt.xlabel(r"Model time $\tau$", fontsize=fsize)
plt.ylabel("Host $h$", fontsize=fsize)
plt.xticks((0, 200, 400, 600, 800, 1000), ("$0$", "$200$", "$400$", "$600$", "$800$", "$1000$"), fontsize=22)
plt.yticks((0.0, 0.2, 0.4, 0.6, 0.8, 1.0), ("$0.0$", "$0.2$", "$0.4$", "$0.6$", "$0.8$", "$1.0$"), fontsize=22)

plt.ylim(0, 1)
for i in xrange(len(gammaList)):
    gamma1 = gammaList[i]
    Dynamics = odeint(System, Initial, TimeList, args=(A, gamma1), rtol=1.e-12,atol=1.e-12)
    #sub1.plot(TimeList, Dynamics[:, 1], color= color_list[i+1], linewidth=4)
    #sub2.plot(TimeList, Dynamics[:, 0], color= color_list[i+1], linewidth=4, label=label_list[i]+"  $\gamma$")
    sub1.plot(TimeList, Dynamics[:, 1], color= color_list[i], linewidth=4)
    sub2.plot(TimeList, Dynamics[:, 0], color= color_list[i], linewidth=4, label=label_list[i]+"  $\gamma$")
        
plt.legend(loc="best", fontsize=19)
plt.xlim((0, Timemax))
fig.savefig('Desktop/Temporal.pdf', transparent=True) # Change this line appropriately when saving file in another directory    
plt.show()    
"""

# Dynamical Trajectories in the phase plane (Ms1 Figure 5 to be runned with previous lines of codes)
"""
LineType = array(["--", "-", "--", "--","--", "--"])

Timemax = 3000
TimeList = linspace(0., Timemax, 20000)

color_list = ((0.482, 0.867, 0.867), (0.259, 0.608, 0.608), (0.145, 0.471, 0.471), (0.031, 0.376, 0.376))  # blue scale
#color_list = ((0.482, 0.867, 0.867), (0.259, 0.608, 0.608), (0.145, 0.471, 0.471), (0.031, 0.376, 0.376))  # blue scale
color_list = ((0.482+0.5, 0.867, 0.867), (0.259+0.5, 0.608-0.1, 0.608-0.1), (0.145+0.5, 0.471-0.1, 0.471), (0.031+0.5, 0.376-0.2, 0.376-0.2))  # red scale
color_list = ((0.898, 0.60, 0.678), (0.857, 0.276, 0.278), (0.333, 0, 0)) # for low, mid and high gamma


plt.ylim(0, 1)
for i in xrange(len(gammaList)):
    gamma1 = gammaList[i]
    FITNESS = FitnessX(UList, XList, gamma1)
    plt.subplot(2,3,i+1+3) #i+1)
    if i+1 == 1:
        plt.text(-1, 0.9, title_list[i], fontsize = 25)
    else:
        plt.text(-0.6, 0.9, title_list[i], fontsize = 25)
    kappa = gamma1**l/(gammaHalfOpt**l + gamma1**l) 
    cost_gam = gamma1**m/(1-gamma1**m)
    FitGrad = -delta*A*UList*(cost_gam)*exp(-delta*UList) - A*(cost_gam)*(1 - exp(-delta*UList)) + kappa*(-XList + 1)*(exp(-UList))
    plt.contour(UList, XList, FITNESS, NullLevels, colors=((0.145, 0.471, 0.471),), linewidths=5)
    plt.contour(UList,XList, FitGrad , NullLevels, colors=((0.482, 0.867, 0.867),), linewidths=5)
    plt.legend()
    plt.xlabel("Trait value $u$", fontsize=fsize)
    if i+1==1:
        plt.ylabel("Host $h$", fontsize = fsize)
    Dynamics = odeint(System, Initial, TimeList, args=(A, gamma1), rtol=1.e-12,atol=1.e-12) 
    #plt.plot(Dynamics[:, 1], Dynamics[:, 0],color=color_list[i+1],linewidth=3, label=label_list[i]+"  $\gamma$")
    plt.plot(Dynamics[:, 1], Dynamics[:, 0],color=color_list[i],linewidth=3, label=label_list[i]+"  $\gamma$")
    plt.plot(Dynamics[len(TimeList)-1, 1], Dynamics[len(TimeList)-1, 0], "o", color=color_list[i], markersize=10)
    plt.xticks((0, 1, 2, 3, 4, 5), ("$0$", "$1$", "$2$", "$3$", "$4$", "$5$"),fontsize=fsize)
    if i+1 ==1:
        plt.yticks((0.0, 0.2, 0.4, 0.6, 0.8, 1.0), ("$0.0$", "$0.2$", "$0.4$", "$0.6$", "$0.8$", "$1.0$"), fontsize=fsize)
    else:
        plt.yticks((0.0, 0.2, 0.4, 0.6, 0.8, 1.0), (" ", " ", " ", " ", " ", " "), fontsize=fsize)
    
      
    plt.legend(loc="upper right", fontsize=19) 

fig.savefig('Desktop/PhasePlane.pdf', transparent=True) # Change this line appropriately when saving file in another directory
              
plt.show()   
"""
# Sensitivity to both A and gamma (Ms1 Figure 7)

AList, gammaList = meshgrid(linspace(0, 0.5, 200), linspace(0, 0.9999999999, 200))
Timemax = 20000
TimeList = linspace(0., Timemax, 15000)

SensitivityH = zeros(AList.shape)
SensitivityU = zeros(AList.shape)

import numpy
for i in xrange(AList.shape[0]):
    for j in xrange(AList.shape[1]):
        print i, j
        A1 = AList[i, j]
        gamma1 = gammaList[i, j]
        subresult = odeint(System, Initial, TimeList, args=(A1, gamma1), rtol=1.e-12,atol=1.e-12)
        SensitivityHpre = mean(subresult[:, 0])
        SensitivityUpre = mean(subresult[:, 1])
        if numpy.isnan(SensitivityUpre) or numpy.isnan(SensitivityHpre): # if we have nan, that mean the numerical solution crashed
            SensitivityH[i, j] = 0.01 # the model crashes because the cost ot symbiosis is too high but we can deduct from the figure that there is a collapse of coral population in those case
            SensitivityU[i, j] = 0.01
        else: 
            SensitivityH[i, j] = SensitivityHpre
            SensitivityU[i, j] = SensitivityUpre

fig=plt.figure(figsize = (25, 10))
plt.subplot(2, 3, 1+3) 
#plt.text(-0.15, 0.9, title_list[0], fontsize = 25)
plt.text(-0.15, 0.9, title_list[0], fontsize = 25)
levels = array([0, 0.5, 1.5, 3, 8])
plt.contour(AList, gammaList, SensitivityU, levels, colors="white", linewidths=3) 
color_list0 = ("#C3B4CB", "#9775AA", "#764B8E", "#3D1255")#)"#260339") # for traits

"""
cs = plt.contourf(AList, gammaList, SensitivityU, levels, colors='none',
                  hatches=['.', "o", "\\", "/"],
                  extend='lower'
                  )# File can only be saved as eps format
""" 
cs = plt.contourf(AList, gammaList, SensitivityU, levels, colors= color_list0,
                  extend='lower'
                  )# File can only be saved as eps format    
artists1, labels_X= cs.legend_elements()

               
labels1 = ["Very low investment", 
                "Low investment", 
                "Mid investment", 
                "High investment"]                 
plt.legend(artists1, labels1, borderpad=0.75, handlelength =3.5, handleheight=1.7, columnspacing=1.6, loc=(-0.05, 1.10), ncol=2, fontsize=15, framealpha = 0)

plt.xticks((0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), ("$0.0$", "$0.1$","$0.2$", "$0.3$","$0.4$", "$0.5$","$0.6$"),fontsize=fsize)
plt.yticks((0.0, 0.2,0.4, 0.6, 0.8, 1.0), ("$0$", "$0.2$","$0.4$", "$0.6$", "$0.8$", "$1.0$"), fontsize=fsize)

plt.xlabel(r"Symbiotic cost parameter $a$", fontsize = fsize)
#plt.xlabel(r"Symbiotic cost parameter $\alpha/\beta G_H  $", fontsize = fsize)

plt.ylabel("Symbiont to host biomass ratio $\gamma$", fontsize = fsize)
#plt.ylabel("Symbiont to host biomass ratio $\Gamma$", fontsize = fsize)
plt.xlim((0, 0.5))

plt.subplot(2, 3, 2+3) 
plt.text(-0.075, 0.9, title_list[1], fontsize = 25)
levels = array([0,  0.02, 0.25, 0.6, 1])

color_list = ((0.482, 0.867, 0.867), (0.259, 0.608, 0.608), (0.145, 0.471, 0.471), (0.031, 0.376, 0.376)) # for corals


plt.contour(AList, gammaList, SensitivityH, levels, colors="white", linewidths=3) 

"""
cs2 = plt.contourf(AList, gammaList, SensitivityH, levels, colors='none',
                  hatches=['.', "o", "\\", "/"],
                  extend='lower'
                  ) # File can only be saved as eps format
""" 
cs2 = plt.contourf(AList, gammaList, SensitivityH, levels, colors=color_list,
                  extend='lower'
                  ) # File can only be saved as eps format    
artists2, labels_X= cs2.legend_elements()

labels2 = ["Collapse", 
                "Low host biomass", 
                "Mid host biomass", 
                "High host biomass"]             
plt.legend(artists2, labels2, borderpad=0.75, handlelength =3.5, handleheight=1.7,columnspacing=1.6, loc=(-0.05, 1.10), ncol=2, fontsize=15, framealpha = 0)

plt.xticks((0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), ("$0.0$", "$0.1$","$0.2$", "$0.3$","$0.4$", "$0.5$","$0.6$"),fontsize=fsize)
plt.yticks((0.0, 0.2,0.4, 0.6, 0.8, 1.0), (" ", " "," ", " ", " ", " "), fontsize=fsize)
plt.xlim((0, 0.5))

#plt.xlabel(r"Cost of symbiosis $a$", fontsize = fsize)
plt.xlabel(r"Symbiotic cost parameter $a$", fontsize = fsize)
#plt.xlabel(r"Symbiotic cost parameter $\alpha/\beta G_H  $", fontsize = fsize)
fig.savefig('Desktop/Sen.pdf', transparent=True) # Change this line appropriately when saving file in another directory

plt.show()


# Shock experiment (Ms1 Figure 8)
"""
def shock(gamma1, gammashock, T, k):
    if T<k*Timemax/3:
        return gamma1
    else:
        return gammashock
def ShockSystem(Initial, T, A, gamma1, gammashock, k):
    dSystem = zeros(len(Initial))
    X = Initial[0]
    U = Initial[1]
    E = 1 - exp(-U)
    GAMMA = shock(gamma1, gammashock, T, k)
    print T
    kappa = GAMMA**l/(gammaHalfOpt**l + GAMMA**l)
    cost_gam = GAMMA**m/(1-GAMMA**m)
    Benefit = kappa*E*(1-X)
    Cost =  (1-exp(-delta*U))*A*U*(cost_gam) + M
    Fitness = (Benefit - Cost)
    dSystem[0] = Fitness*X
    dSystem[1] = SIGMA*Gradient(X, U, A, GAMMA)
    return dSystem


LineType = array(["--", "--", "-", "o", "o","--", "--"])  
Timemax = 1000
TimeList = linspace(0., Timemax, 10000)
fig=plt.figure(figsize = (25, 11))
gammaList = array([0.01, 0.1, 0.3, 0.6, 0.75])  
sub1 = plt.subplot(2, 3, 1)
plt.text(-150, 4.90, title_list[0], fontsize = 25)
#plt.xlabel(r"Scaled time $G_H\,t$", fontsize=fsize)
plt.xlabel(r"Model time $\tau$", fontsize=fsize)
#plt.ylabel(r"Scaled trait value $\beta U$", fontsize=fsize)
plt.ylabel(r"Trait value $U$", fontsize=fsize)
plt.xticks((0, 200, 400, 600, 800, 1000), ("$0$", "$200$", "$400$", "$600$", "$800$", "$1000$"), fontsize=fsize)
plt.yticks((0, 1, 2, 3, 4, 5), ("$0$", "$1$", "$2$", "$3$", "$4$", "$5$"), fontsize=fsize)
plt.ylim((0, 5))
plt.xlim((0, Timemax))  

sub2 = plt.subplot(2, 3, 2)
plt.text(-180, 0.975, title_list[1], fontsize = 25)
plt.xlabel(r"Model time $\tau$", fontsize=fsize)
#plt.ylabel("Scaled coral biomass $H/K_H$", fontsize=fsize)
plt.ylabel("Coral biomass $h$", fontsize=fsize)
plt.xticks((0, 200, 400, 600, 800, 1000), ("$0$", "$200$", "$400$", "$600$", "$800$", "$1000$"), fontsize=fsize)
plt.yticks((0.0, 0.2, 0.4, 0.6, 0.8, 1.0), ("$0.0$", "$0.2$", "$0.4$", "$0.6$", "$0.8$", "$1.0$"), fontsize=fsize)

plt.ylim(0, 1)

#color_list=array([(0.8, 0.8, 0.8),(0.55, 0.55, 0.55),(0, 0, 0), (0.3, 0.3, 0.3), (0, 0, 0)])
#color_list=array(["red","green","black", "blue", "magenta"])
#color_list = array([(0.482, 0.867, 0.867), (0.259, 0.608, 0.608), (0.145, 0.471, 0.471), (0.031, 0.376, 0.376), (0.05, 0.255, 0.255)])
color_list = array([(0.482+0.5, 0.867-0.1, 0.867), (1, 0.255+0.3, 0.255+0.4),(0.259+0.5, 0.608-0.1, 0.608-0.1), (0.145+0.5, 0.471-0.1, 0.471), (0.259+0.2, 0.608-0.3, 0.608-0.3)])
label_list=array(["Severe bleaching", "Moderate bleaching", "No shock", "Moderate outburst", "Severe outburst"])

gamma=0.3
for i in xrange(len(gammaList)):
    gammashock = gammaList[i]
    if i in (0, 1):
        Dynamics = odeint(ShockSystem, Initial, TimeList, args=(A, gamma, gammashock, 1), rtol=1.e-12,atol=1.e-12)
    else:
        Dynamics = odeint(ShockSystem, Initial, TimeList, args=(A, gamma, gammashock, 2), rtol=1.e-12,atol=1.e-12)
    if i in (0, 1):
        sub1.plot(TimeList[TimeList>=1*(Timemax/3)], Dynamics[TimeList>=1*(Timemax/3), 1], linewidth=3, color=color_list[i])#, label="$\gamma$ = %.2f"%gamma1)
        #sub2.plot(TimeList[TimeList>=1*(Timemax/3)], Dynamics[TimeList>=1*(Timemax/3), 0], linewidth=3, color=color_list[i], label=label_list[i]+" ($\gamma_{shock}$ = %d)"%gammashock)
        sub2.plot(TimeList[TimeList>=1*(Timemax/3)], Dynamics[TimeList>=1*(Timemax/3), 0], linewidth=3, color=color_list[i], label=label_list[i])
    elif i in (3, 4):
        sub1.plot(TimeList[TimeList>=2*(Timemax/3)], Dynamics[TimeList>=2*(Timemax/3), 1], linewidth=3, color=color_list[i])#, label="$\gamma$ = %.2f"%gamma1)
        #sub2.plot(TimeList[TimeList>=2*(Timemax/3)], Dynamics[TimeList>=2*(Timemax/3), 0], linewidth=3, color=color_list[i], label=label_list[i]+" ($\gamma_{shock}$ = %d)"%gammashock)
        sub2.plot(TimeList[TimeList>=2*(Timemax/3)], Dynamics[TimeList>=2*(Timemax/3), 0], linewidth=3, color=color_list[i], label=label_list[i])
    else:
        sub1.plot(TimeList, Dynamics[:, 1], linewidth=6, color=color_list[i])#, label="$\gamma$ = %.2f"%gamma1)
        #sub2.plot(TimeList, Dynamics[:, 0], linewidth=6, color=color_list[i], label=label_list[i]+" ($\gamma$ = %d)"%gammashock)
        sub2.plot(TimeList, Dynamics[:, 0], linewidth=6, color=color_list[i], label=label_list[i])
  
plt.legend(loc= "upper right", fontsize=15, framealpha = 0)   
plt.xlim((0, Timemax))  
fig.savefig('Desktop/Shock.pdf', transparent=True) # Change this line appropriately when saving file in another directory
plt.show()
"""
     

    
  

