import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# creating a blank window
# for the animation
fig = plt.figure()

#Better View for UCM
axis = plt.axes(xlim =(-15, 15),ylim =(-15, 15))
fig.set_figwidth(6)
fig.set_figheight(6)

#Better View for Parabolic Path
# axis = plt.axes(xlim =(-50, 50),ylim =(-50, 50))
# fig.set_figwidth(6)
# fig.set_figheight(6)


line, = axis.plot([], [], lw = 2)
lngood, = axis.plot([], [], 'go', markersize=4)

# what will our line dataset
# contain?
def init():
	line.set_data([], [])
	return line,

# initializing empty values
# for x and y co-ordinates
xdata, ydata = [], []

filename=input("Enter the file : ")
f=open(filename,"r")
data=f.readlines()
sz=len(data)
f.close()
# print(data)

def animate(i):
    # print(i)
    myline=data[i].split()
    x=float(myline[0])
    y=float(myline[1])
    # print(x,y)
	# appending values to the previously
	# empty x and y data holders
    xdata.append(x)
    ydata.append(y)
    line.set_data(xdata, ydata)
    lngood.set_data(x,y)
    return line,lngood

# calling the animation function	
anim = animation.FuncAnimation(fig, animate,
							# init_func = init,
							frames = sz,
							interval = 0.1,
							blit = True)
plt.show()
# saves the animation in our desktop
# anim.save('growingCoil.mp4', writer = 'ffmpeg', fps = 30)
print("Completed")