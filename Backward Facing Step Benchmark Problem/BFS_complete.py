# -*- coding: utf-8 -*-

"""

import numpy as np
import matplotlib.pyplot as plt

L1=1;
L2=2;
Lc=1;
L3=10;
dx=dy=0.1;
dt=0.0001;
N=int((L2+L3)/dx)+1;
M=int((L1 +Lc)/dy)+1;
Re=100;
u0=Re/L1;
s_n=np.zeros((N,M));
s_n1=np.zeros((N,M));
w_n=np.zeros((N,M));
w_n1=np.zeros((N,M));
u=np.zeros((N,M));
v=np.zeros((N,M));
it=[];
#Boundary Condition's

#Inlet for u,s_n,s_n1,w_n and w_n1

for j in range (int(Lc/dy),M):
    y=j*dy-Lc;
    u[0][j]=u0*y*(1-(y/L1))/(4*L1);
    s_n[0][j]=s_n1[0][j]=u0*y*y*(0.5-(y/(3*L1)))/(4*L1);
    w_n[0][j]=w_n1[0][j]=-u0*(1-(2*y/L1))/(4*L1);
    
#Upper Wall for u,s_n & s_n1

for i in range (1,N-1):
    u[i][M-1]=0;
    s_n[i][M-1]=s_n1[i][M-1]=u0*L1*L1*(0.5-(1/3))/(4*L1);

#Upper Wall for u,s_n & s_n1 =0 already
   
#Time stepping
   
for t in range (100):
    
    #Upper Wall w_n & w_n1

    for i in range (1,N-1):
        w_n1[i][M-1]=-(2/dy**2)*(s_n[i][M-2]-s_n[i][M-1]);
    
    #Lower inlet Wall w_n1

    for i in range (1,int(L2/dx)+1):
        j=int(Lc/dy);
        w_n1[i][j]=-(2/dy**2)*(s_n[i][j+1]-s_n[i][j]);
    
    #Lower side inlet Wall w_n1

    for j in range (0,int(Lc/dy)+1):
        i=int(L2/dx);
        w_n1[i][j]=-(2/dx**2)*(s_n[i+1][j]-s_n[i][j]);        
        
    #Lower outlet Wall w_n1

    for i in range (int(L2/dx),N-1):
        j=0;
        w_n1[i][j]=-(2/dy**2)*(s_n[i][j+1]-s_n[i][j]);        
    
    # Vorticity equation for upper half
    
    for i in range (1,N-1):
        for j in range (int(Lc/dy)+1,M-1):
            w_n1[i][j]=w_n[i][j]+dt*((w_n[i+1][j]-2*w_n[i][j]+w_n[i-1][j])/(dx*dx)+(w_n[i][j+1]-2*w_n[i][j]+w_n[i][j-1])/(dy*dy)-u[i][j]*(w_n[i+1][j]-w_n[i][j])/(dx)-v[i][j]*(w_n[i][j+1]-w_n[i][j])/(dy));

    # Vorticity equation for Lower half
    
    for i in range (int(L2/dx)+1,N-1):
        for j in range (1,int(Lc/dy)+1):
            w_n1[i][j]=w_n[i][j]+dt*((w_n[i+1][j]-2*w_n[i][j]+w_n[i-1][j])/(dx*dx)+(w_n[i][j+1]-2*w_n[i][j]+w_n[i][j-1])/(dy*dy)-u[i][j]*(w_n[i+1][j]-w_n[i][j])/(dx)-v[i][j]*(w_n[i][j+1]-w_n[i][j])/(dy));
    
    # Vorticity B.C at outlet
    
    for j in range (1,M-1):
        i=N-1;
        w_n1[i][j]=w_n1[i-1][j];
        

    #Add for difference check for stream function
    
    # Solving for stream function
    
    for k in range (30):

 # Stream function equation for upper half
        
        for i in range (1,N-1):
            for j in range (int(Lc/dy)+1,M-1):
                s_n1[i][j]=(w_n1[i][j]*(dx*dx*dy*dy)+(s_n1[i+1][j]+s_n1[i-1][j])*(dy*dy)+(s_n1[i][j+1]+s_n1[i][j-1])*(dx*dx))/(2*(dx*dx+dy*dy));
                
                
# Stream function equation for Lower half
    
        for i in range (int(L2/dx)+1,N-1):
            for j in range (1,int(Lc/dy)+1):
                s_n1[i][j]=(w_n1[i][j]*(dx*dx*dy*dy)+(s_n1[i+1][j]+s_n1[i-1][j])*(dy*dy)+(s_n1[i][j+1]+s_n1[i][j-1])*(dx*dx))/(2*(dx*dx+dy*dy));

# Stream function B.C at outlet
    
        for j in range (1,M):
            i=N-1;
            s_n1[i][j]=s_n1[i-1][j];   
            
        error=0.0
        for i in range(0,N):
            for j in range(0,M):
                error=error+(s_n1[i][j]-s_n[i][j])**2
                s_n[i][j]=s_n1[i][j]          
            error=np.sqrt(error/(((N-int(L2/dx))*(Lc/dy))+(N*(M-Lc/dy))))
            
        if (error)<0.000001:
            break
        it.append(error);    
# solving for u and v
    
    for i in range (1,N):
        for j in range (1,M):
            u[i][j]=(s_n1[i][j]-s_n1[i][j-1])/dy;
            v[i][j]=-(s_n1[i][j]-s_n1[i-1][j])/dx;

#Copying N+1 matrices to N matrices
    
    w_n=w_n1[:];
    s_n=s_n1[:];

length=0;
q=[];
for i in range(int(L2/dx)+1,N):
    j=int(Lc/dy)/2;
    q.append(u[i][j]);
    if (u[i+1][j]-u[i][j]<.01):
        length=i*dx;
        break;

print(length)

plt.plot(q)
    


#plt.plot(it)
#plt.contourf(np.transpose(s_n1),1000);
