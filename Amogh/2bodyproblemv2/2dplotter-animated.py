import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

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

fig=plt.figure()
ax = fig.add_subplot()
ax.set_xlim(-500,500)
ax.set_ylim(-500,500)

line1, =ax.plot([],[],'r',lw=1,animated=True)
line2, =ax.plot([],[],'b',lw=1,animated=True)

ax.set_xlabel('X', fontweight = 'bold', fontsize = 14)
ax.set_ylabel('Y', fontweight = 'bold', fontsize = 14)

def animate(k):
    line1.set_data(x1_coords[:k+1], y1_coords[:k+1])
    line2.set_data(x2_coords[:k+1], y2_coords[:k+1])
    print("Processing frame",k,"/",n_iter)
    
    return [line1, line2]

data_range=[]
for i in range(0,n_iter,100):
    data_range.append(i)
    
anim = animation.FuncAnimation(fig, animate, data_range, interval=100, blit=True, repeat=False)

writergif=animation.PillowWriter(fps=100)

anim.save("simulation.gif", writer=writergif)


plt.show()


