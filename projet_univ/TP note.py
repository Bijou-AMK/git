#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 11:26:33 2023

@author: meili
"""
# w : courant d'adaptation : adapter le degré d'ouverture des canaux ionique (dynamique)
# v : ici v, dans le cours "u" ; correspond au "v" d'un premier modèle simplifié (qui résume les 3 variables n, m, h (on avait commencé par éliminer n puis tout remplacer par minfini)) ; 
# pour le modèle exp on modélise le début de la forme du spike par une focntion exp ; 
# enfin le modèle adex permet de s'adapter au fur et à mesure que se déroule le train

import numpy as np
import math
import matplotlib.pyplot as plt

data ={
'tau_m':20 ,   # constant du temps
'a':0 ,   # la resistance de canal ionique
'tau_w ':100 ,  # constant du temps
'b':5/np . power ( 10 , 12 ) ,
'u_reset ':-60 ,
'u_rest ':-70 ,
'R':500*np . power ( 10 , 9 ) ,  # la resistance de la membrane
'th ':-45 ,    # seuil
'th_rh ':-50 ,  # seuil : tendance a declencher
'DeltaT ':2 ,  # expoentiel
'step ':65/np . power ( 10 , 12 ) # courant injecte
}

T=300       # interval du temps
dt=0.01    # chaque 0.01s on va enregistrer
N=np.floor(T/dt).astype(int) 
I=[data ['step ']*(i*dt>=0) for i in range(N)] # courant injecté (ici toujours 5)

def adex (T , dt ,I , data ):
    N=np . floor ( T/dt ) . astype (int)
    v=np . zeros ( N )
    w=np . zeros ( N )
    
    v[0]= data ['u_rest ']
    DeltaT = data ['DeltaT ']
    th_rh = data ['th_rh ']
    R = data ['R']
    a = data ['a']
    w[0]=0
    
    for i in range ( N-1 ):
        if v[i]> data ['th ']:#si on depasse le seuil
            v[i]=0   # pour voir le spike visuellement
            v[i+1]= data ['u_reset ']
            w[i+1]=w[i]+ data ['b'] # canal plus ouvert
        else :# sinon on resout les EDO
            v[i+1]=v[i]+dt* (((-v[i]+v[0])+ DeltaT * math.exp((v[i]-th_rh)/DeltaT)-R*w[i]+R*I[i]) / data ['tau_m'])
            w[i+1]=w[i]+dt* ((a*(v[i]-v[0])-w[i]) / data ['tau_w '] )# rajouter les formules du modèle (voir cours)
    return [v , w]
   

[v , w] = adex(T,dt,I,data) 


temps = np.linspace(0.0, T, len(v))
plt.plot(temps,v,label='PTM',linewidth=2)
#=================================================================================================
# Qestion1
# Tonic
I=[data ['step ']*(i*dt>=50) for i in range(N)] # courant injecté (ici toujours 5)
data ={
'tau_m':20 ,   # constant du temps
'a':0 ,   # la resistance de canal ionique
'tau_w ':30 ,  # constant du temps
'b':60/np . power ( 10 , 12 ) ,
'u_reset ':-60 ,
'u_rest ':-70 ,
'R':500*np . power ( 10 , 9 ) ,  # la resistance de la membrane
'th ':-45 ,    # seuil
'th_rh ':-50 ,  # seuil : tendance a declencher
'DeltaT ':2 ,  # expentiel
'step ':65/np . power ( 10 , 12 ) # courant injecte
}

[v , w] = adex(T,dt,I,data) 
temps = np.linspace(0.0, T, len(v))
plt.plot(temps,v,linewidth=2)
plt.xlabel("temps (ms)")
plt.ylabel("u(mV)")
plt.title("Tonic")
plt.legend()
##====================================================================
# Bursting
data ={
'tau_m':5 ,   # constant du temps
'a':-0.5/np . power ( 10 , 12 ) ,   # la resistance de canal ionique
'tau_w ':100 ,  # constant du temps
'b':7/np . power ( 10 , 12 ) ,
'u_reset ':-46 ,
'u_rest ':-70 ,
'R':500*np . power ( 10 , 9 ) ,  # la resistance de la membrane
'th ':-45 ,    # seuil
'th_rh ':-50 ,  # seuil : tendance a declencher
'DeltaT ':2 ,  # expentiel
'step ':65/np . power ( 10 , 12 ) # courant injecte
}

[v , w] = adex(T,dt,I,data) 
temps = np.linspace(0.0, T, len(v))
plt.plot(temps,v,linewidth=2)
plt.xlabel("temps (ms)")
plt.ylabel("u(mV)")
plt.title("Bursting")
plt.legend() #plein de spike et des temps de repos 
##=================================================================
## Adapting
data ={
'tau_m':20 ,   # constant du temps
'a':0 ,   # la resistance de canal ionique
'tau_w ':100 ,  # constant du temps
'b':5/np . power ( 10 , 12 ) ,
'u_reset ':-55 ,
'u_rest ':-70 ,
'R':500*np . power ( 10 , 9 ) ,  # la resistance de la membrane
'th ':-45 ,    # seuil
'th_rh ':-50 ,  # seuil : tendance a declencher
'DeltaT ':2 ,  # expentiel
'step ':65/np . power ( 10 , 12 ) # courant injecte
}

[v , w] = adex(T,dt,I,data) 
temps = np.linspace(0.0, T, len(v))
plt.plot(temps,v,linewidth=2)
plt.xlabel("temps (ms)")
plt.ylabel("u(mV)")
plt.title("Adapting")
plt.legend()# ça selargit 
##=================================================================================================
## Qestion 2 
data ={
'tau_m':20 ,   # constant du temps
'a':0 ,   # la resistance de canal ionique
'tau_w ':30 ,  # constant du temps
'b':100/np . power ( 10 , 12 ) ,
'u_reset ':-55 ,
'u_rest ':-70 ,
'R':500*np . power ( 10 , 9 ) ,  # la resistance de la membrane
'th ':-45 ,    # seuil
'th_rh ':-50 ,  # seuil : tendance a declencher
'DeltaT ':2 ,  # expentiel
'step ':65/np . power ( 10 , 12) # courant injecte
}

T=300       # interval du temps
dt=0.01    # chaque 0.01s on va enregistrer
N=np.floor(T/dt).astype(int) 
I=[data ['step ']*(i*dt>=0) for i in range(N)] # courant injecté (ici toujours 5)

def adex (T , dt ,I , data ):
    N=np . floor ( T/dt ) . astype (int)
    v=np . zeros ( N )
    w=np . zeros ( N )
    
    v[0]= data ['u_rest ']
    DeltaT = data ['DeltaT ']
    th_rh = data ['th_rh ']
    R = data ['R']
    a = data ['a']
    w[0]=0
    nb=0
    
    for i in range ( N-1 ):
        if v[i]> data ['th ']:#si on depasse le seuil
            v[i]=0   # pour voir le spike visuellement
            v[i+1]= data ['u_reset ']
            w[i+1]=w[i]+ data ['b'] # canal plus ouvert
        else :# sinon on resout les EDO
            v[i+1]=v[i]+dt* (((-v[i]+v[0])+ DeltaT * math.exp((v[i]-th_rh)/DeltaT)-R*w[i]+R*I[i]) / data ['tau_m'])
            w[i+1]=w[i]+dt* ((a*(v[i]-v[0])-w[i]) / data ['tau_w '] )# rajouter les formules du modèle (voir cours)
        if v[i] == 0 :
            nb += 1
    return [v , w, nb]
   
courant = 0
Id = []
fois = []
while courant < 200/np . power ( 10 , 12 ) :
    I=[courant*(i*dt>=0) for i in range(N)] 
    [v,w,nb]=adex(T,dt,I,data) 
    fois.append(nb)
    Id.append(courant)
    courant += 1/np . power ( 10 , 12 )
plt.plot(Id,fois,label='fonction de gain',linewidth=2)
##=================================================================================================
## Qestion 3
T=50
dt=0.01
N=np.floor(T/dt).astype(int)

Id = []
time = []
courant = 0
td = 0

for td in range(40):
    for courant in np.arange(0,800/np . power ( 10 , 12 ),50/np . power ( 10 , 12 )): 
        
        I=[courant*(i*dt>=10)*(i*dt<=10+td) for i in range(N)] 
        [v,m,nb]=adex(T,dt,I,data) 
        
        if len(time) != 0 and nb == 1 and time[-1]!= td :
                time.append(td)
                Id.append(courant)
        elif nb == 1 and len(time) == 0:
                time.append(td)
                Id.append(courant)
        
plt.plot(time,Id,label='La courbe amplitude-durée',linewidth=2)
plt.xlabel("temps (ms)")
plt.ylabel("amplitude (mV)")
plt.title("La courbe amplitude-durée")        
##=================================================================================================
## Qestion 4
data ={
'tau_m':20 ,   # constant du temps
'a':0 ,   # la resistance de canal ionique
'tau_w ':30 ,  # constant du temps
'b':100/np . power ( 10 , 12 ) ,
'u_reset ':-55 ,
'u_rest ':-70 ,
'R':500*np . power ( 10 , 9 ) ,  # la resistance de la membrane
'th ':-45 ,    # seuil
'th_rh ':-50 ,  # seuil : tendance a declencher
'DeltaT ':2 ,  # expentiel
'step ':65/np . power ( 10 , 12 ) # courant injecte
}
T=300       # interval du temps
dt=0.01    # chaque 0.01s on va enregistrer
N=np.floor(T/dt).astype(int) 
I=[data ['step ']*(i*dt>=0) for i in range(N)] # courant injecté (ici toujours 5)

def adex (T , dt ,I , data ):
    N=np . floor ( T/dt ) . astype (int)
    v=np . zeros ( N )
    w=np . zeros ( N )
    
    v[0]= data ['u_rest ']
    DeltaT = data ['DeltaT ']
    th_rh = data ['th_rh ']
    R = data ['R']
    a = data ['a']
    w[0]=0
    sigma=5
    
    for i in range ( N-1 ):    
        I_alea = sigma*np . random . normal (0 , 1 )/np . power ( 10 , 11 )

        if v[i]> data ['th ']:#si on depasse le seuil
            v[i]=0   # pour voir le spike visuellement
            v[i+1]= data ['u_reset ']
            w[i+1]=w[i]+ data ['b'] # canal plus ouvert
        else :# sinon on resout les EDO
            v[i+1]=v[i]+dt* (((-v[i]+v[0])+ DeltaT * math.exp((v[i]-th_rh)/DeltaT)-R*w[i]+R*(I[i]+I_alea)) / data ['tau_m'])
            w[i+1]=w[i]+dt* ((a*(v[i]-v[0])-w[i]) / data ['tau_w '] )# rajouter les formules du modèle (voir cours)
    return [v , w]
   
for i in range(3):
    [v , w] = adex(T,dt,I,data) 
    temps = np.linspace(0.0, T, len(v))
    plt.plot(temps,v,label='bruit',linewidth=1)
    
#-------------------------------------------------------------------------------------------------
def adex (T , dt ,I , data ):
    N=np . floor ( T/dt ) . astype (int)
    v=np . zeros ( N )
    w=np . zeros ( N )
    
    v[0]= data ['u_rest ']
    DeltaT = data ['DeltaT ']
    th_rh = data ['th_rh ']
    R = data ['R']
    a = data ['a']
    w[0]=0
    sigma=5
    declencher=[]
    
    for i in range ( N-1 ):
        I_alea = sigma*np . random . normal (0 , 1 )/np . power ( 10 , 11 )

        if v[i]> data ['th ']:#si on depasse le seuil
            v[i]=0   # pour voir le spike visuellement
            v[i+1]= data ['u_reset ']
            w[i+1]=w[i]+ data ['b'] # canal plus ouvert
        else :# sinon on resout les EDO
            v[i+1]=v[i]+dt* (((-v[i]+v[0])+ DeltaT * math.exp((v[i]-th_rh)/DeltaT)-R*w[i]+R*(I[i]+I_alea)) / data ['tau_m'])
            w[i+1]=w[i]+dt* ((a*(v[i]-v[0])-w[i]) / data ['tau_w '] )# rajouter les formules du modèle (voir cours)
        if v[i]==0:
            declencher.append(i*dt)
    return [v , w, declencher]    

dixfois = []
for i in range(10):
    [v , w , declencher] = adex(T,dt,I,data) 
    dixfois.append(declencher)
    

y=list(range(1,11))

for xe,ye in zip(dixfois,y):
    plt.scatter([xe] , [ye]*len(xe))