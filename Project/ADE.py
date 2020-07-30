                                                    


import numpy as np
import matplotlib.pyplot as plt
import time

start=time.time();

L1=1;
L2=1;
T1=100;
T2=200;
T3=300;
T4=400;

dx=dy=0.1;
dt=0.00001;

alpha=0.5;

dt_stream = 0.002
epsilon_stream = 1.0e-2
rx=alpha*(dt)/dx**2;
ry=alpha*(dt)/dy**2;

N=int(L1/dx)+1;
M=int(L2/dy)+1;

T_n=np.zeros([N,M]);
T_n1=np.zeros([N,M]);
T_n2=np.zeros([N,M]);

# Boundary Condition's

for j in range (0,M):
    i=0;
    T_n[i][j]=T_n1[i][j]=T_n2[i][j]=T1;
for j in range (0,M):
    i=N-1;
    T_n[i][j]=T_n1[i][j]=T_n2[i][j]=T3;
for i in range (0,N-1):
    j=0;
    T_n[i][j]=T_n1[i][j]=T_n2[i][j]=T2;
for i in range (0,N-1):
    j=M-1;
    T_n[i][j]=T_n1[i][j]=T_n2[i][j]=T4;

 
#Time stepping
er=[];
for t in range (1000):
    for k in range (100):
        
        #step 1
        for i in range (1,N-1):
            for j in range (1,M-1):
                T_n1[i][j]=(rx*T_n1[i-1][j]+(1-ry-rx)*T_n[i][j]+rx*T_n[i+1][j]/
                    +ry*T_n1[i][j-1]+ry*T_n[i][j+1])/(1+rx+ry);
        #Step 2
        for i in range (N-2,0,-1):
            for j in range (M-2,0,-1):
                T_n2[i][j]=(rx*T_n1[i-1][j]+(1-ry-rx)*T_n1[i][j]+rx*T_n2[i+1][j]/
                    +ry*T_n1[i][j-1]+ry*T_n2[i][j+1])/(1+rx+ry);       
        error = 0.0;
        for i in range(0,N):
            for j in range(0,M):
                error = error + (T_n1[i][j] - T_n[i][j])**2
                er.append(error);
                T_n[i][j] = T_n2[i][j]
        error = np.sqrt(error/(M*N))
        print("timestep", t, "error=", error/epsilon_stream, "iteration number", k)
        er.append(error/epsilon_stream);
    if (error/dt_stream)<epsilon_stream:
        break
    
plt.plot(er)
plt.xlabel('Iteration')
plt.ylabel('Error/epsilon_stream')
plt.title('Error per Iteration graph of ADE Scheme')
plt.show()
plt.contourf(T_n)
plt.title('Temperature Distribution using ADE scheme')


end=time.time();
print("Compilation time is: ")
print(end-start)
