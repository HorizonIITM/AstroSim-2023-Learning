import matplotlib.pyplot as plt
import csv
from matplotlib.animation import FuncAnimation

x1,y1=[],[]
x2,y2=[],[]
with open('trial1.csv', mode ='r')as file:
  csvFile = csv.reader(file)
  for lines in csvFile:
    x1.append(float(lines[0]))
    y1.append(float(lines[1]))
    x2.append(float(lines[2]))
    y2.append(float(lines[3]))
plt.style.use('dark_background')
fig = plt.figure()
ax = plt.axes()
ax.set_ylim(-2,2)
ax.set_xlim(-2,2.5)
ax.set_aspect(1)
p1 = plt.Circle((0,0), 2e-2, color="blue")
ax.add_artist(p1)
p2 = plt.Circle((0,0), 2e-2, color="orange")
ax.add_artist(p2)
def anim_func(i):
  ax.plot(x1[:i], y1[:i],color="blue")
  ax.plot(x2[:i], y2[:i],color="orange")
  p1.center = x1[i],y1[i]
  p2.center = x2[i],y2[i]
animation = FuncAnimation(fig, anim_func, frames =  1500,interval = 20)
animation.save("trial1.mp4", dpi=600)