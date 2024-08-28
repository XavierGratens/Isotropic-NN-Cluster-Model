import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

Data = np.loadtxt("C:/Xavier/Xavier2020/Python/MyProgram/NN_Cluster_Model/Quartets/QuartetsXavier.txt", dtype=int)

print(Data.shape[0])
Row=24


print(Data[Row, 1], Data[Row, 2],Data[Row, 3],Data[Row, 4], Data[Row, 5],Data[Row, 6])
print(Data[Row,:])


print(Data[Row, 1], Data[Row, 2],Data[Row, 3],Data[Row, 4], Data[Row, 5],Data[Row, 6])




#print (Data.shape[0])

for i in range(Data.shape[0]):
    for j in range(6):
        if Data[i,j+1]== 2:
            Data[i,j+1] = 1




Gem = np.zeros((4,3))
S1j = np.zeros((3*3,3))
S2j = np.zeros((2*3,3))
S3j = np.zeros((1*3,3))


R=Row

#print(Data[R, 1], Data[R, 2],Data[R, 3],Data[R, 4], Data[R, 5],Data[R, 6])
for i in range(3):
    Gem[0,0]=0
    Gem[0,1]=0
    Gem[0,2]=0
    if Data[R,18+i*3] % 2 ==0:
       Gem[i+1,0]=Data[R,16+i*3]*np.power(3,1/2)/2
       
    else:
       Gem[i+1,0]=Data[R,16+i*3]*np.power(3,1/2)/2 + 1/(np.power(3,1/2)) 
    Gem[i+1,1]=Data[R,17+i*3]-Data[R,16+i*3]*1/2
    Gem[i+1,2]=Data[R,18+i*3]*np.power(2/3,1/2)

print(Gem)

for j in range(3):
        if Data[R,j+1]== 1:
             if Data[R,18+j*3] % 2 ==0:
                 S1j[3*j+1,0]=Data[R,16+j*3]*np.power(3,1/2)/2
             else:
                 S1j[3*j+1,0]=Data[R,16+j*3]*np.power(3,1/2)/2 + 1/(np.power(3,1/2))                      
             S1j[3*j+1,1] = Data[R,17+3*j]-Data[R,16+j*3]*1/2
             S1j[3*j+1,2] = Data[R,18+3*j]*np.power(2/3,1/2)
    
#print(S1j)



def Tube1(RR):
    LR=np.zeros((3,7))    
    for j in range(3):
           if Data[RR,j+1]== 1:
                if Data[RR,18+j*3] % 2 ==0:
                   LR[j,0]=(Data[RR,16+j*3]*np.power(3,1/2)/2)/2
                else:
                   LR[j,0]=(Data[RR,16+j*3]*np.power(3,1/2)/2 + 1/(np.power(3,1/2)))/2                      
                LR[j,1] = (Data[RR,17+3*j]-Data[RR,16+j*3]*1/2)/2
                LR[j,2] = (Data[RR,18+3*j]*np.power(2/3,1/2))/2
                LR[j,3] = np.arccos(LR[j,0]*2)
                LR[j,4] = np.arccos(LR[j,1]*2)
                LR[j,5] = np.arccos(LR[j,2]*2)
                LR[j,6] = np.power(np.power(LR[j,0]*2,2)+np.power(LR[j,1]*2,2)+np.power(LR[j,2]*2,2),1/2)
    return LR

TT1=Tube1(R)


def Tube1A(RR):
    LR=np.zeros((3,6))    
    for j in range(3):
           if Data[RR,j+1]== 1:
                if Data[RR,18+j*3] % 2 ==0:
                   LR[j,0]=(Data[RR,16+j*3]*np.power(3,1/2)/2)/2
                else:
                   LR[j,0]=(Data[RR,16+j*3]*np.power(3,1/2)/2 + 1/(np.power(3,1/2)))/2                      
                LR[j,1] = (Data[RR,17+3*j]-Data[RR,16+j*3]*1/2)/2
                LR[j,2] = (Data[RR,18+3*j]*np.power(2/3,1/2))/2
                LR[j,3] = np.arctan(LR[j,1]/LR[j,0])*180/(np.pi)
                LR[j,4] = np.arccos(LR[j,2]*2)*180/(np.pi)                
                LR[j,5] = np.power(np.power(LR[j,0]*2,2)+np.power(LR[j,1]*2,2)+np.power(LR[j,2]*2,2),1/2)
                if LR[j,0]<0:
                    LR[j,3]=LR[j,3]+180
                if all([LR[j,0]>0,LR[j,1]<0]) :
                   LR[j,3]=LR[j,3]+360    
                    

                
    return LR

TT1=Tube1A(R)


print(TT1)


def Tube2A(RR):
    LR=np.zeros((2,6))
    if Data[RR,18] % 2 ==0:
                   ZX=(Data[RR,16]*np.power(3,1/2)/2)
    else:
                   ZX=(Data[RR,16]*np.power(3,1/2)/2 + 1/(np.power(3,1/2)))                      
    ZY = (Data[RR,17]-Data[RR,16]*1/2)
    ZZ = (Data[RR,18]*np.power(2/3,1/2))
    print(ZX,ZY,ZZ)
    
    for j in range(2):
           if Data[RR,j+1+3]== 1:
                if Data[RR,18+3+j*3] % 2 ==0:
                   LR[j,0]=((Data[RR,16+3+j*3]*np.power(3,1/2)/2)-ZX)/2
                else:
                   LR[j,0]=((Data[RR,16+3+j*3]*np.power(3,1/2)/2 + 1/(np.power(3,1/2)))-ZX)/2                      
                LR[j,1] = ((Data[RR,17+3+3*j]-Data[RR,16+3+j*3]*1/2)-ZY)/2
                LR[j,2] = ((Data[RR,18+3+3*j]*np.power(2/3,1/2))-ZZ)/2

                LR[j,3] = np.arctan(LR[j,1]/LR[j,0])*180/(np.pi)
                LR[j,4] = np.arccos(LR[j,2]*2)*180/(np.pi)                
                LR[j,5] = np.power(np.power(LR[j,0]*2,2)+np.power(LR[j,1]*2,2)+np.power(LR[j,2]*2,2),1/2)
                if LR[j,0]<0:
                    LR[j,3]=LR[j,3]+180
                if all([LR[j,0]>0,LR[j,1]<0]) :
                   LR[j,3]=LR[j,3]+360
                LR[j,0] = LR[j,0]+ZX
                LR[j,1] = LR[j,1]+ZY
                LR[j,2] = LR[j,2]+ZZ
                
                    

                
    return LR

TT2=Tube2A(R)


print(TT2)

def Tube3A(RR):
    LR=np.zeros((1,6))
    if Data[RR,18+3] % 2 ==0:
                   ZX=(Data[RR,16+3]*np.power(3,1/2)/2)
    else:
                   ZX=(Data[RR,16+3]*np.power(3,1/2)/2 + 1/(np.power(3,1/2)))                      
    ZY = (Data[RR,17+3]-Data[RR,16+3]*1/2)
    ZZ = (Data[RR,18+3]*np.power(2/3,1/2))
    print(ZX,ZY,ZZ)
    
    for j in range(1):
           if Data[RR,j+1+5]== 1:
                if Data[RR,18+6+j*3] % 2 ==0:
                   LR[j,0]=((Data[RR,16+6+j*3]*np.power(3,1/2)/2)-ZX)/2
                else:
                   LR[j,0]=((Data[RR,16+6+j*3]*np.power(3,1/2)/2 + 1/(np.power(3,1/2)))-ZX)/2                      
                LR[j,1] = ((Data[RR,17+6+3*j]-Data[RR,16+6+j*3]*1/2)-ZY)/2
                LR[j,2] = ((Data[RR,18+6+3*j]*np.power(2/3,1/2))-ZZ)/2

                LR[j,3] = np.arctan(LR[j,1]/LR[j,0])*180/(np.pi)
                LR[j,4] = np.arccos(LR[j,2]*2)*180/(np.pi)                
                LR[j,5] = np.power(np.power(LR[j,0]*2,2)+np.power(LR[j,1]*2,2)+np.power(LR[j,2]*2,2),1/2)
                if LR[j,0]<0:
                    LR[j,3]=LR[j,3]+180
                if all([LR[j,0]>0,LR[j,1]<0]) :
                   LR[j,3]=LR[j,3]+360
                LR[j,0] = LR[j,0]+ZX
                LR[j,1] = LR[j,1]+ZY
                LR[j,2] = LR[j,2]+ZZ
                
                    

                
    return LR

TT3=Tube3A(R)


print(TT3)




for k in range(2*3):
      if Data[R,18] % 2 ==0:
       S2j[k,0]=Data[R,16]*np.power(3,1/2)/2
      else:
       S2j[k,0]=Data[R,16]*np.power(3,1/2)/2 + 1/(np.power(3,1/2))                      
      S2j[k,1] = Data[R,17]-Data[R,16]*1/2
      S2j[k,2] = Data[R,18]*np.power(2/3,1/2)


for j in range(2):
        if Data[R,j+4]== 1:
             if Data[R,21+j*3] % 2 ==0:
                 S2j[3*j+1,0]=Data[R,19+j*3]*np.power(3,1/2)/2
             else:
                 S2j[3*j+1,0]=Data[R,19+j*3]*np.power(3,1/2)/2 + 1/(np.power(3,1/2))                      
             S2j[3*j+1,1] = Data[R,20+3*j]-Data[R,19+j*3]*1/2
             S2j[3*j+1,2] = Data[R,21+3*j]*np.power(2/3,1/2)

#print(S2j)

for k in range(1*3):
      if Data[R,21] % 2 ==0:
       S3j[k,0]=Data[R,19]*np.power(3,1/2)/2
      else:
       S3j[k,0]=Data[R,19]*np.power(3,1/2)/2 + 1/(np.power(3,1/2))                      
      S3j[k,1] = Data[R,20]-Data[R,19]*1/2
      S3j[k,2] = Data[R,21]*np.power(2/3,1/2)

for j in range(1):
        if Data[R,j+6]== 1:
             if Data[R,24+j*3] % 2 ==0:
                 S3j[3*j+1,0]=Data[R,22+j*3]*np.power(3,1/2)/2
             else:
                 S3j[3*j+1,0]=Data[R,22+j*3]*np.power(3,1/2)/2 + 1/(np.power(3,1/2))                      
             S3j[3*j+1,1] = Data[R,23+3*j]-Data[R,22+j*3]*1/2
             S3j[3*j+1,2] = Data[R,24+3*j]*np.power(2/3,1/2)


#print(S3j)






fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

ax.set_xlim(-2,2)
ax.set_ylim(-2,2) 
ax.set_zlim(-2,2) 
ax.set_aspect('equal')


    
ax.scatter(Gem[:,0], Gem[:,1], Gem[:,2])
ax.plot(S1j[:,0], S1j[:,1], S1j[:,2])
ax.plot(S2j[:,0], S2j[:,1], S2j[:,2])
ax.plot(S3j[:,0], S3j[:,1], S3j[:,2])
plt.show()



np.savetxt("C:/Xavier/Xavier2020/Python/MyProgram/NN_Cluster_Model/Quartets/QuartetsXavier.txt", Data,fmt="%d", delimiter=' ' )
