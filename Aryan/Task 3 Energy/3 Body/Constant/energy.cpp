#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

const double G = 1; // Gravitational constant
// const double h = 0.01;    // Time step
// const double T = 20000.0; // Total time

// Class for body
class Body {
public:
  double mass;
  // double energy;
  vector<double> position;
  vector<double> momentum;

  // constructor
  Body(double m, double x, double y, double z, double px, double py,
       double pz /* ,double E */) {
    mass = m;
    position = {x, y, z};
    momentum = {px, py, pz};
  }
};

// Class for Runge Kutta Solver
class Runge_Kutta {
private:
  // PosX derivative
  vector<double> x_dash(vector<Body> B, int n) {
    vector<double> dx(n, 0); // Initialising all the elements to 0
    for (int i = 0; i < n; i++) {
      dx[i] = B[i].momentum[0] /
              B[i].mass; // x_dash = dx/dt = vx = px/m ... // momenetum[0] = px
    }
    return dx;
  }

  // PosY derivative
  vector<double> y_dash(vector<Body> B, int n) {
    vector<double> dy(n, 0); // Initialising all the elements to 0
    for (int i = 0; i < n; i++) {
      dy[i] = B[i].momentum[1] /
              B[i].mass; // y_dash = dy/dt = vy = py/m ... // momenetum[1] = py
    }
    return dy;
  }

  // PosZ derivative
  vector<double> z_dash(vector<Body> B, int n) {
    vector<double> dz(n, 0); // Initialising all the elements to 0
    for (int i = 0; i < n; i++) {
      dz[i] = B[i].momentum[2] /
              B[i].mass; // z_dash = dz/dt = vz = pz/m ... // momenetum[2] = pz
    }
    return dz;
  }

  // MomentumX derivative
  vector<double> px_dash(vector<Body> B, int n) {
    vector<double> dpx(n, 0); // Initialising all the elements to 0

    // Selecting pairs to calculate force between them
    for (int i = 0; i < n - 1; i++) {   // First body of pair
      for (int j = i + 1; j < n; j++) { // Second body of pair

        // Calculating distance between 2 bodies
        double rx =
            B[j].position[0] -
            B[i].position[0]; // Position of 2nd body - 1st body in X-Drxn
        double ry =
            B[j].position[1] -
            B[i].position[1]; // Position of 2nd body - 1st body in Y-Drxn
        double rz =
            B[j].position[2] -
            B[i].position[2]; // Position of 2nd body - 1st body in Z-Drxn
        double r =
            sqrt(rx * rx + ry * ry + rz * rz); // Distance between 2 bodies

        // Force F = G*m1*m2/r^2, acceleration of m1 ax = G*m2*rx/r^3
        // dpx/dt = m*dvx/dt = m*ax = G*m1*m2*rx/r^3 ... negative for other body
        // due to equal and opposite force
        double force = (G * B[i].mass * B[j].mass) / (r * r * r);

        dpx[i] += force * rx;
        dpx[j] -= force * rx;
      }
    }

    return dpx;
  }

  // MomentumY derivative
  vector<double> py_dash(vector<Body> B, int n) {
    vector<double> dpy(n, 0); // Initialising all the elements to 0

    // Selecting pairs to calculate force between them
    for (int i = 0; i < n - 1; i++) {   // First body of pair
      for (int j = i + 1; j < n; j++) { // Second body of pair

        // Calculating distance between 2 bodies
        double rx =
            B[j].position[0] -
            B[i].position[0]; // Position of 2nd body - 1st body in X-Drxn
        double ry =
            B[j].position[1] -
            B[i].position[1]; // Position of 2nd body - 1st body in Y-Drxn
        double rz =
            B[j].position[2] -
            B[i].position[2]; // Position of 2nd body - 1st body in Z-Drxn
        double r =
            sqrt(rx * rx + ry * ry + rz * rz); // Distance between 2 bodies

        // Force F = G*m1*m2/r^2, acceleration of m1 ay = G*m2*ry/r^3
        // dpy/dt = m*dvy/dt = m*ay = G*m1*m2*ry/r^3 ... negative for other body
        // due to equal and opposite force
        double force = (G * B[i].mass * B[j].mass) / (r * r * r);

        dpy[i] += force * ry;
        dpy[j] -= force * ry;
      }
    }

    return dpy;
  }

  // MomentumZ derivative
  vector<double> pz_dash(vector<Body> B, int n) {
    vector<double> dpz(n, 0); // Initialising all the elements to 0

    // Selecting pairs to calculate force between them
    for (int i = 0; i < n - 1; i++) {   // First body of pair
      for (int j = i + 1; j < n; j++) { // Second body of pair

        // Calculating distance between 2 bodies
        double rx =
            B[j].position[0] -
            B[i].position[0]; // Position of 2nd body - 1st body in X-Drxn
        double ry =
            B[j].position[1] -
            B[i].position[1]; // Position of 2nd body - 1st body in Y-Drxn
        double rz =
            B[j].position[2] -
            B[i].position[2]; // Position of 2nd body - 1st body in Z-Drxn
        double r =
            sqrt(rx * rx + ry * ry + rz * rz); // Distance between 2 bodies

        // Force F = G*m1*m2/r^2, acceleration of m1 az = G*m2*rz/r^3
        // dpy/dt = m*dvz/dt = m*az = G*m1*m2*rz/r^3 ... negative for other body
        // due to equal and opposite force
        double force = (G * B[i].mass * B[j].mass) / (r * r * r);

        dpz[i] += force * rz;
        dpz[j] -= force * rz;
      }
    }

    return dpz;
  }

  // Finding Runge-Kutta Constants
  vector<Body> k1(vector<Body> B, int n, double h) {
    vector<double> dx1 = x_dash(B, n);
    vector<double> dy1 = y_dash(B, n);
    vector<double> dz1 = z_dash(B, n);
    vector<double> dpx1 = px_dash(B, n);
    vector<double> dpy1 = py_dash(B, n);
    vector<double> dpz1 = pz_dash(B, n);

    vector<Body> k1_body = B;

    // Updating Values
    for (int i = 0; i < n; i++) {
      k1_body[i].position[0] += h * dx1[i] / 2;  // x2
      k1_body[i].position[1] += h * dy1[i] / 2;  // y2
      k1_body[i].position[2] += h * dz1[i] / 2;  // z2
      k1_body[i].momentum[0] += h * dpx1[i] / 2; // px2
      k1_body[i].momentum[1] += h * dpy1[i] / 2; // py2
      k1_body[i].momentum[2] += h * dpz1[i] / 2; // py2
    }
    return k1_body;
  }

  vector<Body> k2(vector<Body> B, int n, double h, vector<Body> k1_body) {
    vector<double> dx2 = x_dash(k1_body, n);
    vector<double> dy2 = y_dash(k1_body, n);
    vector<double> dz2 = z_dash(k1_body, n);
    vector<double> dpx2 = px_dash(k1_body, n);
    vector<double> dpy2 = py_dash(k1_body, n);
    vector<double> dpz2 = pz_dash(k1_body, n);

    vector<Body> k2_body = B;

    // Updating Values
    for (int i = 0; i < n; i++) {
      k2_body[i].position[0] += h * dx2[i] / 2;  // x3
      k2_body[i].position[1] += h * dy2[i] / 2;  // y3
      k2_body[i].position[2] += h * dz2[i] / 2;  // z3
      k2_body[i].momentum[0] += h * dpx2[i] / 2; // px3
      k2_body[i].momentum[1] += h * dpy2[i] / 2; // py3
      k2_body[i].momentum[2] += h * dpz2[i] / 2; // py3
    }
    return k2_body;
  }

  vector<Body> k3(vector<Body> B, int n, double h, vector<Body> k2_body) {
    vector<double> dx3 = x_dash(k2_body, n);
    vector<double> dy3 = y_dash(k2_body, n);
    vector<double> dz3 = z_dash(k2_body, n);
    vector<double> dpx3 = px_dash(k2_body, n);
    vector<double> dpy3 = py_dash(k2_body, n);
    vector<double> dpz3 = pz_dash(k2_body, n);

    vector<Body> k3_body = B;

    // Updating Values
    for (int i = 0; i < n; i++) {
      k3_body[i].position[0] += h * dx3[i];  // x4
      k3_body[i].position[1] += h * dy3[i];  // y4
      k3_body[i].position[2] += h * dz3[i];  // z4
      k3_body[i].momentum[0] += h * dpx3[i]; // px4
      k3_body[i].momentum[1] += h * dpy3[i]; // py4
      k3_body[i].momentum[2] += h * dpz3[i]; // py4
    }
    return k3_body;
  }

  vector<Body> k4(vector<Body> B, int n, double h, vector<Body> k3_body) {
    vector<double> dx4 = x_dash(k3_body, n);
    vector<double> dy4 = y_dash(k3_body, n);
    vector<double> dz4 = z_dash(k3_body, n);
    vector<double> dpx4 = px_dash(k3_body, n);
    vector<double> dpy4 = py_dash(k3_body, n);
    vector<double> dpz4 = pz_dash(k3_body, n);

    vector<Body> k4_body = B;

    // Updating Values
    for (int i = 0; i < n; i++) {
      k4_body[i].position[0] += h * dx4[i];
      k4_body[i].position[1] += h * dy4[i];
      k4_body[i].position[2] += h * dz4[i];
      k4_body[i].momentum[0] += h * dpx4[i];
      k4_body[i].momentum[1] += h * dpy4[i];
      k4_body[i].momentum[2] += h * dpz4[i];
    }
    return k4_body;
  }

public:


  vector<Body> rk4(vector<Body> B, int n, double h) {
    vector<Body> k1_body = k1(B, n, h);
    vector<Body> k2_body = k2(B, n, h, k1_body);
    vector<Body> k3_body = k3(B, n, h, k1_body);
    vector<Body> k4_body = k4(B, n, h, k1_body);

    for (int i = 0; i < n; i++) {
      B[i].position[0] +=
          (2 * k1_body[i].position[0] + 4 * k2_body[i].position[0] +
           2 * k3_body[i].position[0] + k4_body[i].position[0] -
           9 * B[i].position[0]) /
          6;
      B[i].position[1] +=
          (2 * k1_body[i].position[1] + 4 * k2_body[i].position[1] +
           2 * k3_body[i].position[1] + k4_body[i].position[1] -
           9 * B[i].position[1]) /
          6;
      B[i].position[2] +=
          (2 * k1_body[i].position[2] + 4 * k2_body[i].position[2] +
           2 * k3_body[i].position[2] + k4_body[i].position[2] -
           9 * B[i].position[2]) /
          6;
      B[i].momentum[0] +=
          (2 * k1_body[i].momentum[0] + 4 * k2_body[i].momentum[0] +
           2 * k3_body[i].momentum[0] + k4_body[i].momentum[0] -
           9 * B[i].momentum[0]) /
          6;
      B[i].momentum[1] +=
          (2 * k1_body[i].momentum[1] + 4 * k2_body[i].momentum[1] +
           2 * k3_body[i].momentum[1] + k4_body[i].momentum[1] -
           9 * B[i].momentum[1]) /
          6;
      B[i].momentum[2] +=
          (2 * k1_body[i].momentum[2] + 4 * k2_body[i].momentum[2] +
           2 * k3_body[i].momentum[2] + k4_body[i].momentum[2] -
           9 * B[i].momentum[2]) /
          6;
    }
    return B;
  }

  double totalEnergy(vector<Body> B, int n) {
    double totalEnergy = 0.0;

    for (int i = 0; i < n; i++)

    {
      double kineticEnergy = 0.5 *
                             (B[i].momentum[0] * B[i].momentum[0] +
                              B[i].momentum[1] * B[i].momentum[1] +
                              B[i].momentum[2] * B[i].momentum[2]) /
                             B[i].mass;

      double potentialEnergy = 0.0; 

      for (int j = i + 1; j < n; j++) {
        double rx = B[j].position[0] - B[i].position[0];
        double ry = B[j].position[1] - B[i].position[1];
        double rz = B[j].position[2] - B[i].position[2];
        double r = sqrt(rx * rx + ry * ry + rz * rz);

        potentialEnergy -= (G * B[i].mass * B[j].mass) / r;
      }

      totalEnergy += kineticEnergy + potentialEnergy;
    }
    return totalEnergy;
  }
};

int main() {
  // Initial Conditions
  vector<Body> B = {Body(100, 0, 0, 0, 0, 0, 0), Body(1, 16, 0, 0, 0, 2.5, 0),
                    Body(1, 25, 0, 0, 0, 2.449, 0)};
  double T = 10;    // Total time
  double h = 0.00001; // Time step

  Runge_Kutta gravity_solver;
  ofstream file1, energyFile;

  file1.open("positions.txt");
  energyFile.open("energy_3.txt");

  for (double t = 0; t <= T; t += h) {
    B = gravity_solver.rk4(B, 3, h);
    double totalEnergy = gravity_solver.totalEnergy(B, 3);

    file1 << B[0].position[0] << "," << B[0].position[1] << ","
          << B[0].position[2] << "," << B[1].position[0] << ","
          << B[1].position[1] << "," << B[1].position[2] << ","
          << B[2].position[0] << "," << B[2].position[1] << ","
          << B[2].position[2] << '\n';

    energyFile << t << " " << totalEnergy << '\n';
  }

  file1.close();
  energyFile.close();
  cout << "Completed" << endl;
  return 0;
}
