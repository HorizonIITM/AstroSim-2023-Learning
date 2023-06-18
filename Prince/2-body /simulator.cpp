#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

const double G = 6.67430e-11;

struct Body {
    double mass;
    double x, y;
    double vx, vy;
};

void updateAcceleration(Body& heavyBody, Body& lightBody) {
    double dx = lightBody.x - heavyBody.x;
    double dy = lightBody.y - heavyBody.y;
    double distanceSquared = dx * dx + dy * dy;
    double invDistance = 1.0 / std::sqrt(distanceSquared);
    double force = G * heavyBody.mass * invDistance * invDistance;
    double ax = force * dx * invDistance / lightBody.mass;
    double ay = force * dy * invDistance / lightBody.mass;
    lightBody.vx -= ax;
    lightBody.vy -= ay;
}

void updatePosition(Body& body, double timeStep) {
    body.x += body.vx * timeStep;
    body.y += body.vy * timeStep;
}

int main() {
    // Creating two bodies
    Body heavyBody, lightBody;

    // Properties of the bodies
    heavyBody.mass = 1.0e6;
    heavyBody.x = 0.0;
    heavyBody.y = 0.0;
    heavyBody.vx = 0.0;
    heavyBody.vy = 0.0;

    lightBody.mass = 1.0;
    lightBody.x = 100.0;
    lightBody.y = 100.0;
    lightBody.vx = 0.0;
    lightBody.vy = -10.0;

    // Simulation parameters
    double t_start = 0.0;
    double t_end = 10.0;
    double timeStep = 0.01;

    // Output file
    std::ofstream outputFile("output.csv");
    if (!outputFile) {
        std::cout << "Error opening file." << std::endl;
        return 1;
    }

    // Simulating the motion
    for (double t = t_start; t <= t_end; t += timeStep) {
        // Update the acceleration of the light body due to the heavy body
        updateAcceleration(heavyBody, lightBody);

        // Updating the position of the light body
        updatePosition(lightBody, timeStep);

        // Output
        outputFile << lightBody.x << "," << lightBody.y << std::endl;
    }

    outputFile.close();
    std::cout << "Simulation completed." << std::endl;
    return 0;
}
