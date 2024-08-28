import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

Data = np.loadtxt("Triplets.txt")

print(Data)
print(Data[0, 1], Data[0, 2],Data[0, 3])

#print (Data.shape[0])

for i in range(Data.shape[0]):
    for j in range(3):
        if Data[i,j+1]== 2:
            Data[i,j+1] = 1

print(Data)


Gem = np.zeros((3,3))
S1j = np.zeros((2*3,3))
S2j = np.zeros((1*3,3))



R=6

for i in range(2):
    Gem[0,0]=0
    Gem[0,1]=0
    Gem[0,2]=0
    if Data[R,12+i*3] % 2 ==0:
       Gem[i+1,0]=Data[R,10+i*3]*np.power(3,1/2)/2
       
    else:
       Gem[i+1,0]=Data[R,10+i*3]*np.power(3,1/2)/2 + 1/(np.power(3,1/2)) 
    Gem[i+1,1]=Data[R,11+i*3]-Data[R,10+i*3]*1/2
    Gem[i+1,2]=Data[R,12+i*3]*np.power(2/3,1/2)

print(Gem)


for j in range(2):
        if Data[R,j+1]== 1:
             if Data[R,12+j*3] % 2 ==0:
                 S1j[3*j+1,0]=Data[R,10+j*3]*np.power(3,1/2)/2
             else:
                 S1j[3*j+1,0]=Data[R,10+j*3]*np.power(3,1/2)/2 + 1/(np.power(3,1/2))                      
             S1j[3*j+1,1] = Data[R,11+3*j]-Data[R,10+j*3]*1/2
             S1j[3*j+1,2] = Data[R,12+3*j]*np.power(2/3,1/2)

print(S1j)

for k in range(1*3):
      if Data[R,12] % 2 ==0:
       S2j[k,0]=Data[R,10]*np.power(3,1/2)/2
      else:
       S2j[k,0]=Data[R,10]*np.power(3,1/2)/2 + 1/(np.power(3,1/2))                      
      S2j[k,1] = Data[R,11]-Data[R,10]*1/2
      S2j[k,2] = Data[R,12]*np.power(2/3,1/2)

for j in range(1):
        if Data[R,j+3]== 1:
             if Data[R,15+j*3] % 2 ==0:
                 S2j[3*j+1,0]=Data[R,13+j*3]*np.power(3,1/2)/2
             else:
                 S2j[3*j+1,0]=Data[R,13+j*3]*np.power(3,1/2)/2 + 1/(np.power(3,1/2))                      
             S2j[3*j+1,1] = Data[R,14+3*j]-Data[R,13+j*3]*1/2
             S2j[3*j+1,2] = Data[R,15+3*j]*np.power(2/3,1/2)

print(S2j)



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')


    
ax.scatter(Gem[:,0], Gem[:,1], Gem[:,2])
ax.plot(S1j[:,0], S1j[:,1], S1j[:,2])
ax.plot(S2j[:,0], S2j[:,1], S2j[:,2])
#ax.plot(S3j[:,0], S3j[:,1], S3j[:,2])
plt.show()
