import matplotlib.pyplot as plt
import csv
from matplotlib.animation import FuncAnimation

x,y=[],[]

with open('output.csv', mode ='r')as file:
  csvFile = csv.reader(file)
  for lines in csvFile:
    x.append(float(lines[0]))
    y.append(float(lines[1]))
xp=x[::1000]
yp=y[::1000]

plt.style.use('dark_background')
fig = plt.figure()
ax = plt.axes()
ax.set_ylim(-1e12,1e12)
ax.set_xlim(-1e12,1e12)
ax.set_aspect(1)

earth = plt.Circle((0,0), 5e9, color="blue")
ax.add_artist(earth)
sun = plt.Circle((0,0), 1e10, color="orange")
ax.add_artist(sun)
def anim_func(i):
  ax.plot(xp[:i], yp[:i],color="blue")
  earth.center = xp[i],yp[i]

animation = FuncAnimation(fig, anim_func, frames =  306,interval = 20)
animation.save("hyperbolic.mp4", dpi=600)
