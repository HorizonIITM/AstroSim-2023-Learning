#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

const double G = 6.674e-11;

const double M = 1e30; // Heavy point mass
const double m = 1; // Light point mass

double Ax(double Rx, double Ry){
  double R = sqrt(Rx*Rx + Ry*Ry);
  double acc = -(G*M)/(R*R*R);
  return acc*Rx;
}

  double Ay(double Rx, double Ry){
  double R = sqrt(Rx*Rx + Ry*Ry);
  double acc = -(G*M)/(R*R*R);
  return acc*Ry;
  
}

void RungeKutta(double& Rx, double& Ry, double& Vx, double& Vy, double h){

  double K1x, K1y, K2x, K2y, K3x, K3y, K4x, K4y; // RK values for r
  double L1x, L1y, L2x, L2y, L3x, L3y, L4x, L4y; // RK values for v

  

  K1x = Vx * h;
  K1y = Vy * h;
  L1x = Ax(Rx,Ry) * h;
  L1y = Ay(Rx,Ry) * h;

  double R2x = Rx + 0.5*K1x;
  double R2y = Ry + 0.5*K1y;
  double V2x = Vx + 0.5*L1x;
  double V2y = Vy + 0.5*L1y;


  K2x = V2x * h;
  K2y = V2y * h;
  L2x = Ax(R2x,R2y) * h;
  L2y = Ay(R2x,R2y) * h;

  double R3x = Rx + 0.5*K2x;
  double R3y = Ry + 0.5*K2y;
  double V3x = Vx + 0.5*L2x;
  double V3y = Vy + 0.5*L2y;


  K3x = V3x * h;
  K3y = V3y * h;
  L3x = Ax(R3x,R3y) * h;
  L3y = Ay(R3x,R3y) * h;

  double R4x = Rx + K3x;
  double R4y = Ry + K3y;
  double V4x = Vx + L3x;
  double V4y = Vy + L3y;


  K4x = V4x * h;
  K4y = V4y * h;
  L4x = Ax(R4x,R4y) * h;
  L4y = Ay(R4x,R4y) * h;

  Rx += (K1x + K2x + K3x + K4x)/6;
  Ry += (K1y + K2y + K3y + K4y)/6;
  Vx += (L1x + L2x + L3x + L4x)/6;
  Vy += (L1y + L2y + L3y + L4y)/6;
}

int main(){
  double Rx, Ry, Vx, Vy;

  // Initial Conditions
  Rx = 1e10;
  Ry = 0;
  Vx = -1e5;
  Vy = 1e5;

  double h = 100; //time step
  double totaltime = 365*24*3600; //total time of simulation

  ofstream outfile("x-y coordinates.txt"); // to store the points

  for (double t = 0; t <= totaltime; t+= h)
    {
      outfile << Rx <<" "<< Ry << endl;

      RungeKutta(Rx, Ry, Vx, Vy, h);
    }

  outfile.close();

  return 0;
  
  
}