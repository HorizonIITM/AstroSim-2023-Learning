

G = 6.67430e-11;    
dt = 0.1;
m1=2e30;
m2=2e20;        
t_end=input("Enter the end time for the simulation:");
x1=0;
y1=0;
vx1=0;
vy1=0;
x2=input("Enter the X cordinate of the lighter object:");
y2=input("Enter the Y cordinate of the lighter object:");
vx2=input("Enter the velocity along X of the lighter object:");
vy2=input("Enter the velocity along Y of the lighter object:");

X1(1)=x1;
Y1(1)=y1;
X2(1)=x2;
Y2(1)=y2;
VX2(1)=vx2;
VY2(1)=vy2;
t = 0;
n=t_end/0.01;

for i=1:n
    r_squared = (X2(i) - x1)^2 + (Y2(i) - y1)^2;
    r = sqrt(r_squared);
    F = G * (m1 * m2) / r_squared;

    Ax2 = F * (x1 - X2(i)) / (m2 * r);
    Ay2 = F * (y1 - Y2(i)) / (m2 * r);
    
    k1x = Ax2 * dt;
    k1y = Ay2 * dt;
    k2x = (Ax2 + k1x/2) * dt;
    k2y = (Ay2 + k1y/2) * dt;
    k3x = (Ax2 + k2x/2) * dt;
    k3y = (Ay2 + k2y/2) * dt;
    k4x = (Ax2 + k3x) * dt;
    k4y = (Ay2 + k3y) * dt;
    VX2(i+1) = VX2(i) + (k1x + 2*k2x + 2*k3x + k4x)/6;
    VY2(i+1) = VY2(i) + (k1y + 2*k2y + 2*k3y + k4y)/6;
    
    k1x = VX2(i+1) * dt;
    k1y = VY2(i+1) * dt;
    k2x = (VX2(i+1) + k1x/2) * dt;
    k2y = (VY2(i+1) + k1y/2) * dt;
    k3x = (VX2(i+1) + k2x/2) * dt;
    k3y = (VY2(i+1) + k2y/2) * dt;
    k4x = (VX2(i+1) + k3x) * dt;
    k4y = (VY2(i+1) + k3y) * dt;
    X2(i+1) = X2(i) + (k1x + 2*k2x + 2*k3x + k4x)/6;
    Y2(i+1) = Y2(i) + (k1y + 2*k2y + 2*k3y + k4y)/6;

    t = t + dt;
end

plot(X2, Y2)
