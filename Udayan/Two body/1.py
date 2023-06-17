from matplotlib import pyplot as plt
import numpy as np
import math as ma

G = 6.67e-11

# the position of heavy mass is at origin (0,0)
div = 10
dt = 86400/div
t_end = 86400 * 365 * 20
#inital masses and the velocity of both bodies
masses = [1.989e30, 1898.19e24]
velocities = np.array([[0, 0], [0, 18e3]])
#the inital position of the of heavy body is 1 and lighter weight is 2
x1 = [0]
y1 = [0]
x2 = [778.570e9]
y2 = [0]
m1 = masses[0]
m2 = masses[1]
i =0
t =0
while (t <= t_end):
    r = ma.sqrt((x2[i]-x1[i])**2 + (y2[i]-y1[i])**2)
    Fx = (G*m1*m2)*(x2[i]-x1[i])/r**3
    Fy = (G*m1*m2)*(y2[i]-y1[i])/r**3
    velocities[0][0] += -1*Fx*dt/masses[0]
    velocities[0][1] += -1*Fy*dt/masses[0]
    velocities[1][0] += -1*Fx*dt/masses[1]
    velocities[1][1] += -1*Fy*dt/masses[1]

    x1.append(x1[i]+ velocities[0][0]*dt)
    y1.append(y1[i]+ velocities[0][1]*dt)
    x2.append(x2[i]+ velocities[1][0]*dt)
    y2.append(y2[i]+ velocities[1][1]*dt)
    #print(y1[i])
    t = dt + t
    i = i + 1

fig, (ax1, ax2) = plt.subplots(2, 1)

# Plot the first dataset on the first subplot
ax1.plot(x1, y1, color='blue', label='Heavy Mass')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_title('Plot 1')
ax1.legend()

# Plot the second dataset on the second subplot
ax2.plot(x2, y2, color='red', label='Light Mass')
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_title('Plot 2')
ax2.legend()

# Adjust the spacing between subplots
plt.tight_layout()

# Display the figure
plt.show()

scale_factor = 1e-9  # Adjust the scale factor as needed

plt.scatter(np.array(x1[1:]) * scale_factor, np.array(y1[1:]) * scale_factor, color='blue', label='Body 1', s=2)
plt.scatter([x1[0] * scale_factor], [y1[0] * scale_factor], color='blue', marker='o', label='Initial Position 1')
plt.scatter(np.array(x2) * scale_factor, np.array(y2) * scale_factor, color='red', label='Body 2', s=2)
plt.scatter([x2[0] * scale_factor], [y2[0] * scale_factor], color='red', marker='o', label='Initial Position 2')
plt.xlabel('X position (in units of 1e9 meters)')
plt.ylabel('Y position (in units of 1e9 meters)')
plt.title('Positions of Two Bodies')
plt.axis('equal')
plt.legend()
plt.show()