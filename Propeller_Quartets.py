import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA



"Constants"

"Bohr magneton in cgs"
muB_cgs=9.2741E-21
"Planck Constant in SI"
hp=6.6262E-34

"Boltzman´s constant in cgs"
kB_cgs=1.3807E-16
" Avogadro´s constant "
Na = 6.022E23 


Spin=3/2





"Landé Factors"
g=2





"Making matrix"

def MSZ1(S):
    Mat=np.zeros(shape=(int(np.power(S*2+1,4)),int(np.power(S*2+1,4))))
    for j in range(int(np.power(S*2+1,4))):
            Mat[j,j] = S-int(j/((2*S+1)*(2*S+1)*(2*S+1)))
    return Mat    
    
def MSZ2(S):
    Mat=np.zeros(shape=(int(np.power(S*2+1,4)),int(np.power(S*2+1,4))))
    for j in range(int(np.power(S*2+1,4))):
            Mat[j,j] = S-int(j/((2*S+1)*(2*S+1)))+(2*S+1)*int(j/((2*S+1)*(2*S+1)*(2*S+1)))
    return Mat    

def MSZ3(S):
    Mat=np.zeros(shape=(int(np.power(S*2+1,4)),int(np.power(S*2+1,4))))
    for j in range(int(np.power(S*2+1,4))):
            Mat[j,j] = S-int(j/((2*S+1)))+(2*S+1)*int(j/((2*S+1)*(2*S+1)))
    return Mat

def MSZ4(S):
    Mat=np.zeros(shape=(int(np.power(S*2+1,4)),int(np.power(S*2+1,4))))
    for j in range(int(np.power(S*2+1,4))):
            Mat[j,j] = S-j+(2*S+1)*int(j/((2*S+1)))
    return Mat

def MSPlus1(S):
    Mat=np.zeros(shape=(int(np.power(S*2+1,4)),int(np.power(S*2+1,4))))
    for i in range(int(np.power(S*2+1,4))):
        for j in range(int(np.power(S*2+1,4))):
            if i+(2*S+1)*(2*S+1)*(2*S+1)==j:
               M = S-int(j/((2*S+1)*(2*S+1)*(2*S+1))) 
               Mat[i,j] = np.power((S*(S+1)-(M)*(M+1)),1/2)
    return Mat    

def MSMinus1(S):
    Mat=np.zeros(shape=(int(np.power(S*2+1,4)),int(np.power(S*2+1,4))))
    for i in range(int(np.power(S*2+1,4))):
        for j in range(int(np.power(S*2+1,4))):
            if i-(2*S+1)*(2*S+1)*(2*S+1)==j:
               M = S-int(j/((2*S+1)*(2*S+1)*(2*S+1))) 
               Mat[i,j] = np.power((S*(S+1)-(M)*(M-1)),1/2)
    return Mat

def MSPlus2(S):
    Mat=np.zeros(shape=(int(np.power(S*2+1,4)),int(np.power(S*2+1,4))))
    for i in range(int(np.power(S*2+1,4))):
        for j in range(int(np.power(S*2+1,4))):
            if i+(2*S+1)*(2*S+1)==j:
               M = S-int(j/((2*S+1)*(2*S+1)))+(2*S+1)*int(j/((2*S+1)*(2*S+1)*(2*S+1)))
               Mat[i,j] = np.power((S*(S+1)-M*(M+1)),1/2)
    return Mat

def MSMinus2(S):
    Mat=np.zeros(shape=(int(np.power(S*2+1,4)),int(np.power(S*2+1,4))))
    for i in range(int(np.power(S*2+1,4))):
        for j in range(int(np.power(S*2+1,4))):
            if i-(2*S+1)*(2*S+1)==j:
               M = S-int(j/((2*S+1)*(2*S+1)))+(2*S+1)*int(j/((2*S+1)*(2*S+1)*(2*S+1)))
               Mat[i,j] = np.power((S*(S+1)-M*(M-1)),1/2)
    return Mat 

def MSPlus3(S):
    Mat=np.zeros(shape=(int(np.power(S*2+1,4)),int(np.power(S*2+1,4))))
    for i in range(int(np.power(S*2+1,4))):
        for j in range(int(np.power(S*2+1,4))):
            if i+(2*S+1)==j:
               M =  S-int(j/((2*S+1)))+(2*S+1)*int(j/((2*S+1)*(2*S+1)))
               Mat[i,j] = np.power((S*(S+1)-M*(M+1)),1/2)
    return Mat

def MSMinus3(S):
    Mat=np.zeros(shape=(int(np.power(S*2+1,4)),int(np.power(S*2+1,4))))
    for i in range(int(np.power(S*2+1,4))):
        for j in range(int(np.power(S*2+1,4))):
            if i-(2*S+1)==j:
               M =  S-int(j/((2*S+1)))+(2*S+1)*int(j/((2*S+1)*(2*S+1)))
               Mat[i,j] = np.power((S*(S+1)-M*(M-1)),1/2)
    return Mat  

def MSPlus4(S):
    Mat=np.zeros(shape=(int(np.power(S*2+1,4)),int(np.power(S*2+1,4))))
    for i in range(int(np.power(S*2+1,4))):
        for j in range(int(np.power(S*2+1,4))):
            if i==j-1:
               M =  S-j+(2*S+1)*int(j/((2*S+1)))
               Mat[i,j] = np.power((S*(S+1)-M*(M+1)),1/2)
    return Mat

def MSMinus4(S):
    Mat=np.zeros(shape=(int(np.power(S*2+1,4)),int(np.power(S*2+1,4))))
    for i in range(int(np.power(S*2+1,4))):
        for j in range(int(np.power(S*2+1,4))):
            if i==j+1:
               M =  S-j+(2*S+1)*int(j/((2*S+1)))
               Mat[i,j] = np.power((S*(S+1)-M*(M-1)),1/2)
    return Mat

SZ1=MSZ1(Spin)
SZ2=MSZ2(Spin)
SZ3=MSZ3(Spin)
SZ4=MSZ4(Spin)
SP1=MSPlus1(Spin)
SP2=MSPlus2(Spin)
SP3=MSPlus3(Spin)
SP4=MSPlus4(Spin)

SM1=MSMinus1(Spin)
SM2=MSMinus2(Spin)
SM3=MSMinus3(Spin)
SM4=MSMinus4(Spin)

SX1=(SP1+SM1)/2
SX2=(SP2+SM2)/2
SX3=(SP3+SM3)/2
SX4=(SP4+SM4)/2

SY1=(SP1-SM1)/2*-1j
SY2=(SP2-SM2)/2*-1j
SY3=(SP3-SM3)/2*-1j
SY4=(SP4-SM4)/2*-1j



EXPropeller=np.dot(SX1,SX2)+np.dot(SY1,SY2)+np.dot(SZ1,SZ2)+np.dot(SX2,SX3)+np.dot(SY2,SY3)+np.dot(SZ2,SZ3)+np.dot(SX2,SX4)+np.dot(SY2,SY4)+np.dot(SZ2,SZ4)



"Zeeman energy as function of the Landé factor, theta and phi angle "
def Zeeman(g,theta,phi):
    ZZ=g*(np.sin(np.pi*theta/180)*np.sin(np.pi*phi/180)*(SX1+SX2+SX3+SX4)+np.sin(np.pi*theta/180)*np.cos(np.pi*phi/180)*(SY1+SY2+SY3+SY4)+np.cos(np.pi*theta/180)*(SZ1+SZ2+SZ3+SZ4))
    return ZZ




"Energy Levels in Kelvin, Field in kOe"
def EvsH_cgs(g,theta,phi,Jex,Field,NP):
   v0 = np.zeros(int((2*Spin+1)*(2*Spin+1)*(2*Spin+1)*(2*Spin+1)+1))
   for i in range(NP+2):
       v=LA.eigvals((muB_cgs/kB_cgs) * Zeeman(g,theta,phi) * i /NP*Field*1000 + 2*Jex*EXPropeller)
       v=v.real
       v=np.sort(v)
       R=np.hstack((i /NP*Field,v))
       v0=np.vstack((v0,R))
   return v0[1:]    

def M_Quartets(g,theta,phi,Jex,Field,NP,Temp):
    E = EvsH_cgs(g,theta,phi,Jex,Field,NP)
    Zero=[0,0]
    Em=E.transpose()
    Em1=Em[1:]
    Em2=Em1.transpose()
    H= Em[0].transpose()
    S = np.zeros(shape=(Em2.shape[0],2))
    mag = np.zeros(shape=(S.shape[0]-1,2))
    MFF=np.zeros(shape=(NP+1,2))
    for i in range(Em2.shape[0]):
        S[i,1] = np.log(np.sum(((np.exp(-(Em2[i]-(np.min(Em2[i])))/Temp)))))-np.min(Em2[i])/Temp
        S[i,0] = H[i]
    for i in range(S.shape[0]-1):
        mag[i,1] = ((S[i+1,1])-(S[i,1]))/((H[1]-H[0])*1000)*Temp/4*kB_cgs/muB_cgs/(2*Spin)
        mag[i,0] = H[i] +(H[1]-H[0])/2
    M = np.vstack((Zero,mag))
    "Interpolation"
    x_i=np.arange(NP+1)
    x_i =x_i/NP*Field
    y_i=np.interp(x_i,M[:,0],M[:,1])
    MFF[:,0] = x_i
    MFF[:,1] = y_i
    return MFF
    
def N_D(Mat):
    Result=np.zeros(shape=(Mat.shape[0]-1,2))
    RI=np.zeros(shape=(Mat.shape[0]-1,2))    
    for i in range (Mat.shape[0]-1):
        Result[i,1]=(Mat[i+1,1]-Mat[i,1])/(Mat[i+1,0]-Mat[i,0])
        Result[i,0]= Mat[i,0] +(Mat[i+1,0]-Mat[i,0])/2
    "Interpolation"
    x_i=np.arange(Mat.shape[0]-1)
    x_i =x_i/(Mat.shape[0]-1)*np.max(Mat[:,0])
    y_i=np.interp(x_i,Result[:,0],Result[:,1])
    RI[:,0] = x_i
    RI[:,1] = y_i    
    return RI    


Jex=0.5
Temp1 = 0.02


Energy=EvsH_cgs(g,0,0,Jex,120,1000)
M1=M_Quartets(g,0,0,Jex,120,1000,Temp1)
Deriv1=N_D(M1)

Energy[:,0]=Energy[:,0]*muB_cgs/kB_cgs*1000*g/(2*Jex)

M1[:,0]=M1[:,0]*muB_cgs/kB_cgs*1000*g/(2*Jex)
Deriv1[:,0]=Deriv1[:,0]*muB_cgs/kB_cgs*1000*g/(2*Jex)
Deriv1[:,1]=Deriv1[:,1]/(muB_cgs/kB_cgs*1000*g/(2*Jex))

fig1=plt.figure(1)
plt.xlabel('$g\mu_{B}H/2 k_B|J_{ex}|$')
plt.ylabel('$Energy/2|J_{ex}|$')

plt.text(0, np.max(Energy[:,Energy.shape[1]-1]*2/3) , 'AF Propeller Quartets, Spin = %s ' % Spin , fontsize=10)

"Plot of the energy levels as function of the field"
for i in range (int((2*Spin+1)*(2*Spin+1)*(2*Spin+1)*(2*Spin+1))):
    plt.plot(Energy[:,0],Energy[:,i+1])
    


fig2=plt.figure(2)
plt.xlabel('$H_ r = g\mu_{B}H/2 k_B|J_{ex}|$')
plt.ylabel('$m = M/M_0$')

plt.text(np.max(M1[:,0])*1/3, np.max(M1[:,1])/8, 'AF Propeller Quartets, Spin = %s ' % Spin , fontsize=10)


"Plot of the energy levels as function of the field"
plt.plot(M1[:,0],M1[:,1], label='$T / |J_{ex}|$ = %s ' % float(Temp1/Jex))
plt.legend()


fig3=plt.figure(3)
plt.xlabel('$H_r = g\mu_{B}H/2 k_B|J_{ex}|$')
plt.ylabel('$dm/dH_r$')

plt.plot(Deriv1[:,0],Deriv1[:,1], label='$T / |J_{ex}|$ = %s ' % float(Temp1/Jex))
plt.legend()


plt.show()

