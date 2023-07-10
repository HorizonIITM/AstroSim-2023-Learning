#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

const double G = 1; // Gravitational constant
const double h = 0.01; // Time step
const double T = 20000.0; // Total time

// Class for Gravitational system
class GravitationalSystem {
public:
    double M; // Heavy mass
    double m; // Light mass
    double x, y; // Positions of the light mass
    double px, py; // Momenta of the light mass

    // Constructor
    GravitationalSystem(double heavyMass, double lightMass, double positionX, double positionY, double momentumX, double momentumY)
        : M(heavyMass), m(lightMass), x(positionX), y(positionY), px(momentumX), py(momentumY) {}

    // Function for defining differential equations
    void derivativesCalculation(double& x_dash, double& y_dash, double& px_dash, double& py_dash) {
        double r = sqrt(x * x + y * y);
        double force = -G * M * m / (r * r * r);

        x_dash = px / m; // Derivative of x position = dx/dt = vx = px / m
        y_dash = py / m; // Derivative of y position = dy/dt = vy = py / m
        px_dash = force * x; // Derivative of x momentum = dpx/dt = m * dvx/dt = m * ax = -GMmx / r^3
        py_dash = force * y; // Derivative of y momentum = dpy/dt = m * dvy/dt = m * ay = -GMmy / r^3
    }

    // Function for Runge-Kutta Integration
    void RungeKutta(double& x, double t, double& dx) {
        double k1 = h * dx;
        double k2 = h * (dx + 0.5 * k1);
        double k3 = h * (dx + 0.5 * k2);
        double k4 = h * (dx + k3);

        x += (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
        t += h;
    }
};

int main() {
    // Initial conditions
    double heavyMass = 100.0; // Heavy point mass
    double lightMass = 1.0; // Light point mass
    double x0 = 20.0; // Initial x position
    double y0 = 0.0; // Initial y position
    double px0 = 0.0; // Initial x momentum
    double py0 = sqrt(double(1.5) * G * heavyMass / x0); // Initial y momentum

    GravitationalSystem system(heavyMass, lightMass, x0, y0, px0, py0);
    double x_dash, y_dash, px_dash, py_dash;

    ofstream outfile("positions.txt");
    if (!outfile) {
        cerr << "Error opening output file." << endl;
        return 1;
    }

    for (double t = 0; t <= T; t += h) {
        outfile << system.x << " " << system.y << endl;

        system.derivativesCalculation(x_dash, y_dash, px_dash, py_dash);

        system.RungeKutta(system.x, t, x_dash);
        system.RungeKutta(system.y, t, y_dash);
        system.RungeKutta(system.px, t, px_dash);
        system.RungeKutta(system.py, t, py_dash);
    }

    outfile.close();
    return 0;
}
