import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

Data = np.loadtxt("Quartets.txt", dtype=int)

print(Data.shape[0])
Row=0
print(Data[Row, 1], Data[Row, 2],Data[Row, 3],Data[Row, 4], Data[Row, 5],Data[Row, 6])
print(Data[Row,:])
Data[1, 5] = 0
Data[1, 6] = 1

Data[5, 5] = 0
Data[5, 6] = 1

Data[7, 5] = 0
Data[7, 6] = 2

Data[9, 5] = 0
Data[9, 6] = 2

Data[18, 5] = 0
Data[18, 6] = 1

Data[45, 5] = 0
Data[45, 6] = 1

Data[46, 5] = 0
Data[46, 6] = 1

Data[47, 5] = 0
Data[47, 6] = 1

Data[48, 5] = 0
Data[48, 6] = 1

Data[49, 5] = 0
Data[49, 6] = 1

Data[50, 5] = 0
Data[50, 6] = 1

Data[51, 5] = 0
Data[51, 6] = 1

Data[56, 3] = 0
Data[56, 6] = 1

print(Data[Row, 1], Data[Row, 2],Data[Row, 3],Data[Row, 4], Data[Row, 5],Data[Row, 6])


np.savetxt("QuartetsXavier.txt", Data,fmt="%d", delimiter=' ' )


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
    
print(S1j)

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

print(S2j)

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


print(S3j)






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



#np.savetxt("C:/Xavier/Xavier2020/Python/MyProgram/NN_Cluster_Model/Quartets/QuartetsXavier.txt", Data,fmt="%d", delimiter=' ' )
