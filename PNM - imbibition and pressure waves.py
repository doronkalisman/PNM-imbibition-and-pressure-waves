# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 20:38:54 2023

@author: doronka
"""


import random
import numpy as np
import scipy.special
import matplotlib.pyplot as plt
import scipy.optimize as opt
import array as arr
from scipy.stats import rankdata
import matplotlib.colors as mcolors


#van-Genuchten parameters 
th_r = 0.03
th_s = 0.39
alfa = 1.1#[1/m] , must be a float
n_ = 3.0 # must be a float
m_ = 1-(1/n_)
l_vg=0.5
Ks = 8e-5 #m/s

#water properties
sigma = 0.072#surface tension [N/m]
vis = 1.002*10**-3#dynamic viscosity [N s/m2]
ro = 1000#density [kg/m3]
g = 9.8#gravity [m/s2]
nu = 1.0*10**-6#kinematic viscosity [m^2/s]
###contact angle
alpha=0#np.pi/3
###    

    
#capillary bundle properties
discrete=25.#Discretezation of capillary groups
dth = (th_s-th_r)/discrete #delta theta increments
const = dth*sigma**2/(2*vis*ro*g)#constant term without tau (tortuosity)

#constructing the nubmer of segments of bundels (M)
M1 = (th_s-th_r)/dth
M2 = np.floor(M1)#number of theta steps
if M2 == M1:
    M = int(M2)
else:
    M = int(M1)+1
loc_v = np.arange(0,M,1)
loc1_v = np.arange(1,M,1)

#determining VG-Mualem retention curve
Se_v = np.nan*np.ones(loc_v.shape)
Se_v[0] = 0.999
Se_v[1:M] = ((th_s-loc1_v*dth)-th_r)/(th_s-th_r)
theta_v = th_r+Se_v*(th_s-th_r)
h_v = ((1/Se_v**(1/m_)-1)**(1/n_))/alfa #calculating h from VG
cap_r_v = 2*sigma/(ro*g*h_v)
#number radius [per m^2]

#Horizontal size of network
junctions=300#must be even and much larger than discretezation (at least oreder of magnitude)

se0=0.55# Dictates the amount of water at each time step as fraction of layer size, requires calibration with CF data

distance=18.# [cm] Longitudinal distance of wetting if completely saturated (addressing porosity)
rounds=int(distance*1.5)


dis_unit=0.005 # Minimum distance unit (m)
dis_rat=int(0.01/dis_unit) # Ratio between minimum distance unit and 1 cm (0.01 m)   

water_mass=int(np.round(se0*junctions/2,0))# Volume of water entering the system at each time step proportional to horizontal dimensions of the system (assuming mainly 1D water propagatin in the vertical direction)

must_reps=int(np.round(dis_rat*distance/(water_mass/(junctions/2.)),0))+1 # The amount of repetitions required to 

layer=must_reps+1# how many layers of the mesh to assure all repetitions have room

num_cap_v=[] #number of capillaries from each group per unit area
#with each group maintaining the same surface area (smaller capillaries=>more capillaries) 
a=0
while a<len(cap_r_v):
    num_cap_v =np.append(num_cap_v,dth/(np.pi*cap_r_v[a]**2))
    a=a+1
    

bigc=2 #number of different realizations of domain

#### Result vectors
bigm=np.zeros((rounds*dis_rat,bigc))
big1=np.zeros((rounds*dis_rat,bigc))

bigtr=np.zeros((rounds*dis_rat,bigc))

bigP1=np.zeros((rounds*dis_rat,bigc))
bigP2=np.zeros((rounds*dis_rat,bigc))
bigP3=np.zeros((rounds*dis_rat,bigc))
bigp=np.zeros((rounds*dis_rat,bigc))

#####
#### FLOW DIRECTION
teta=0   #np.pi 
    #0 for downward flow, pi for upward flow, pi/2 for horizontal
    
    


## Loop of domain realizations begins.
AA=0
while AA<bigc:
    
    Se=[] # Setting the water content entity
    
    SEV=[se0]
    
    
    track=[]
    TR=[]
    

    

##################constructing the geometry of a specific realization
    junc=np.zeros((layer+1,int(junctions/2)))
    cap=np.zeros((layer+1,int(junctions/2)))
    cap_h=np.zeros((layer+1,int(junctions/2)))
    
    cap_rand=np.zeros((layer,int(junctions/2+1)))
    cap_rand_h=np.zeros((layer,int(junctions/2+1)))
    
    cap_ratio=np.round(num_cap_v/sum(num_cap_v)*(junctions)*(layer+1),0)
    
    
    cap_list=[]
    
    a=0
    while a<len(cap_ratio):
        b=0
        while b<cap_ratio[a]:
            cap_list=np.append(cap_list,cap_r_v[a])
            b=b+1
        a=a+1
    
    rand_list=[]
    rand_list=random.sample(list(cap_list),len(cap_list))
    
    cap_rand0=[]    
    cap_rand0=random.sample(list(rand_list[0:int(np.round(len(cap_list)))]),len(cap_list))
    
    a=0
    while a<layer:
        cap_rand[a,:]=cap_rand0[int(a*junctions/2):int((a+1)*junctions/2+1)]
     
        a=a+1
    
    cap_rand0_h=[]    
    cap_rand0_h=random.sample(rand_list[int(np.round(len(cap_list)/2)):],int(len(cap_list)/2))
    
    a=0
    while a<layer:
        cap_rand_h[a,:]=cap_rand0_h[int(a*junctions/2):int((a+1)*junctions/2+1)]
    
        a=a+1
    
#################################################  
    
    ## loop of a single domain realization
    A=0
    while A<must_reps:#   
        fill=[]
                
        
        
        fit=3.3e-4 # [-]  attenuation factor of pressure wave, soil specific, not including radii distribution (in next section)
        ## this parameter symbols the mechanical properties that govern the wave velocity. In some soils, this can be a function of water content.
        ## 1e-4 to 1e-3 is the order of magnitude 
        
        lamda_v=[]
        
        
        ### pressure wave at magnitude at the boundary, can vary over the course of the simulation
        if A<39:
            P_rel_mf=[0e5]# [bar]
        else:
            P_rel_mf=[0e5]
            
        ### The pressure propagation. Determining the location/saturation dependency of the attenuation coefficient and afterwards the normalized pressure magnitude of the wave over distance. 
        
        a=0
        while a<A+1:
            lamda_v=np.ones(len(cap_r_v))
            if a==0:
                lamda_v=fit/cap_r_v[int(np.round((1-se0)*len(cap_r_v))):]
            else:
                lamda_v=fit/cap_r_v[int(np.round((1-SEV[a])*len(cap_r_v))):]
            
            P_rel_per_rad_mf=np.ones(len(cap_r_v))
            P_rel_per_rad_mf=P_rel_mf[a]*np.exp(-lamda_v*dis_unit)
            
                      
            if len(P_rel_per_rad_mf)!= 0:
                P_rel_mf=np.append(P_rel_mf, sum(P_rel_per_rad_mf)/len(P_rel_per_rad_mf))
            else:
                P_rel_mf=np.append(P_rel_mf,0)
                
            a=a+1
        
        
        ### Summing the pressure at each capillary connection, vertical (cap) or horizontal (cap_h)
        
        a=0
        while a<A+1:
            
            for n in range(int(junctions/2)):
                
                cap[a,n]=2.*sigma*np.cos(alpha)/(2.*sigma*np.cos(alpha)/cap_rand[a,n]+(dis_unit*a*np.cos(teta)*ro*g)+P_rel_mf[a])
                cap_h[a,n]=2.*sigma*np.cos(alpha)/(2.*sigma*np.cos(alpha)/cap_rand_h[a,n]+(dis_unit*a*np.cos(teta)*ro*g)+P_rel_mf[a])
            
            a=a+1
            
            
        ### Identifying the water-air interface. This is where the magic happens. Different loops are needed at the beginning of the simulation.
        
        front1=[]
        front2=[]
        
        front1_h=[]
        front2_h=[]
        
        if A==0:
            front1=np.zeros(int(junctions/2))
            front2=np.arange(junctions/2)
            
        if A==1:
            a=0
            while a<junctions/2:
               
                if junc[0,a]!=1:
                    front1.append(0)
                    front2.append(a)
                
                else:
                        
                    if junc[0,a-1]==0:
                    
                            front1_h.append(0)
                            front2_h.append(a-1)
                    
                    if a==junctions/2-1:
                        if junc[0,0]==0:
                            front1_h.append(0)
                            front2_h.append(0)
                            
                      
                    if a!=junctions/2-1:
                        
                        if junc[0,a+1]==0:
                            
                            front1_h.append(0)
                            front2_h.append(a+1)
                    
                    front1.append(1)
                    front2.append(a)
                    
                a=a+1
                            
        if A>1:
                
            a=0
            while a<junctions/2:
               
                if junc[0,a]!=1:
                    front1.append(0)
                    front2.append(a)
                
                b=0
                while b<A:
                    if junc[b,a]==1:
                        if junc[b,a-1]==0:
                    
                            front1_h.append(b)
                            front2_h.append(a-1)
                        
                        if a==junctions/2-1:
                            if junc[b,0]==0:
                                front1_h.append(b)
                                front2_h.append(0)
                            
                      
                        if a!=junctions/2-1:
                        
                            if junc[b,a+1]==0:
                            
                                front1_h.append(b)
                                front2_h.append(a+1)
                                
                        if junc[b+1,a]==0:
                            
                            front1.append(b+1)
                            front2.append(a)
                        
                        if b!=0:
                        
                            if junc[b-1,a]==0:
                            
                                front1.append(b-1)
                                front2.append(a)
                                
                        b=b+1
                    
                    else:
                        b=b+1
                a=a+1
                
        cap_front=[]
        
        limit=1
          
        b=0
        while b<len(front1):
            if cap[int(front1[b]),int(front2[b])]<limit:
                cap_front.append(cap[int(front1[b]),int(front2[b])])
            
            b=b+1
            
        if A!=0: 
            b=0
            while b<len(front1_h):
                if cap_h[int(front1_h[b]),int(front2_h[b])]<limit:
                    cap_front.append(cap_h[int(front1_h[b]),int(front2_h[b])])
                
                b=b+1
            
        ## ranking the capillaries at the interface accorging to pressure magnitude
        
        cap_front0=[]
        cap_front1=[]
        e=[]
        
        cap_front0=(rankdata(cap_front,method='ordinal') - 1).astype(int)
        cap_front1=(rankdata(cap_front,method='min')).astype(int)
      
        e=np.asarray(np.where(cap_front < cap_front[int(np.asarray(np.where(cap_front0 == water_mass)))]))
    
        
        ## Filling the ranked capillaries with water until the determined water volume of the iteration is reached. It's SHOWTIME 
        
        fill=[]
        burst=[]
        burst1=[]
        
        for i, j in enumerate(cap_front1):
                if j < len(e[0]):
                    fill=np.append(fill,i)
                elif len(e[0]) <= j <= water_mass:
                    burst=np.append(burst,i)
        
        burst_list=list(burst)            
        burst1=random.sample(burst_list, water_mass-len(e[0]))
        fill=np.append(fill,burst1)


        ## Updating the domain matrix on water content
        
        for n in range(len(fill)):
            if fill[n]<len(front1):
                junc[int(front1[int(fill[n])]),int(front2[int(fill[n])])]=1.
            else:
                junc[int(front1_h[int(fill[n])-len(front1)]),front2_h[int(fill[n]-len(front1))]]=1.
           
       
         ## Calculation of relative water content per layer 
        
        SEV=np.zeros(layer+1)
        for n in range(len(junc)):
            SEV[n]=float(np.count_nonzero(junc[n,:] == 1))/(junctions/2)
         
        
        ## Time series of water distribution
        Se.append(SEV)
        
        ###################################################
      
               
        A=A+1

   
    
    ###Summing results of iteration
    
    for i in range(rounds*dis_rat):
        bigm[i,AA]=Se[len(Se)-1][i]
        
    
    for i in range(rounds*dis_rat):
        big1[i,AA]=np.average(Se[len(Se)-1][i*2*dis_rat:(i+1)*2*dis_rat-1])
        
    AA=AA+1
    
###props for results summary and display
biga=np.zeros(rounds*dis_rat)
big2=np.zeros(rounds)
bigst2=np.zeros(rounds)
sta=np.zeros(rounds*dis_rat)
upbound=np.zeros(rounds*dis_rat)
downbound=np.zeros(rounds*dis_rat)
        
for i in range(rounds*dis_rat):
    biga[i]=np.average(bigm[i])
    sta[i]=np.std(bigm[i], ddof=1)


    upbound[i]=biga[i]+sta[i]
    downbound[i]=biga[i]-sta[i]
    

for i in range(rounds):
    big2[i]=np.average(big1[i])
    bigst2[i]=np.std(big1[i], ddof=1)
    

    
#%%

####Display##########
    
plt.figure(1,figsize=(10,10)) ### Water content as effective saturation (SE),spatial average and model SD

plt.fill_betweenx(-np.arange(rounds*dis_rat)*dis_unit,upbound,downbound,alpha=0.5, edgecolor='pink', facecolor='pink')
plt.plot(biga,-np.arange(rounds*dis_rat)*dis_unit,color='deepskyblue')
plt.errorbar(big2,-(np.arange(rounds)*2+1)*dis_unit*dis_rat,xerr=bigst2,yerr=0,fmt='o',color='red',markersize=10)

plt.xlabel('Se (-)',fontsize=20)
plt.ylabel('Distance (m)',fontsize=20)
 
plt.xticks(size=14)
plt.yticks(size=14)


plt.ylim((rounds+1)*dis_unit*(-1)*dis_rat-dis_unit,dis_unit)

########


plt.figure(3,figsize=(24,6)) ### Network visualation (final stage) 
figloc=np.arange(len(junc))*dis_unit*(-1)

for i in np.arange(len(junc[0])):
    plt.plot([np.arange(len(junc[0]))[i],np.arange(len(junc[0]))[i]],[(layer)*dis_unit*(-1)*dis_rat,0],color='lightgray')
    
    
for i in np.arange(len(junc)):
    plt.plot([-1,junctions/2],[figloc[i],figloc[i]],color='lightgray')
    plt.plot(np.arange(len(junc[0])),junc[i]*figloc[i],"o",color='blue')

plt.xlim(-1,junctions/2)    
plt.ylim((rounds+1)*dis_unit*(-1)*dis_rat-dis_unit,dis_unit)
    
   

    

