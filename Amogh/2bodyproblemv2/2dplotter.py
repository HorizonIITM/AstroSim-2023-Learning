import numpy as np
from matplotlib import pyplot as plt

n_iter=50000
n_bodies=2

x1_coords=[]
y1_coords=[]

x2_coords=[]
y2_coords=[]

loc_data=open("data_pos.txt","r")
 
for j in range(0,n_iter):
    input_list=loc_data.readline().strip().split()
    x1_coords.append(float(input_list[0]))
    y1_coords.append(float(input_list[1]))


loc_data.readline()

for j in range(0,n_iter):
    input_list=loc_data.readline().strip().split()
    x2_coords.append(float(input_list[0]))
    y2_coords.append(float(input_list[1]))

    
loc_data.close()

plt.figure(1)
plt.plot(x1_coords,y1_coords)
plt.show()

plt.figure(2)
plt.plot(x2_coords,y2_coords)
plt.show()

