data = readtable('draw.txt');
x = data.x;
y = data.y;

plot(0, 0, 'bo');
hold on;
plot(x(1), y(1), 'ro');
plot(x, y, 'r');
plot(x(end), y(end), 'ro');
hold off;

xlabel('x');
ylabel('y');
title('Plot');
legend('Start Point', 'End Point', 'Path');
