% Read the data from the output file
data = csvread('output.csv');
light_x = data(:, 1);
light_y = data(:, 2);

% Plot the heavy object
heavy_x = 0;
heavy_y = 0;
plot(heavy_x, heavy_y, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
hold on

% Create the plot for the lighter body
light_plot = plot(light_x(1), light_y(1), 'bo', 'MarkerSize', 5, 'LineWidth', 1);

% Set the axis limits
axis equal
xlim([-2000, 2000])
ylim([-2000, 2000])

% Animation loop
for i = 2:numel(light_x)
    % Update the position of the lighter body in the plot
    light_plot.XData = light_x(i);
    light_plot.YData = light_y(i);
    
    % Pause for a short duration to visualize the animation
    pause(0.01)
end
