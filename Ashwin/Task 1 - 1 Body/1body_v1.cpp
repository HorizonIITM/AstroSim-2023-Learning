#include <bits/stdc++.h>
#include <vector>
using namespace std;
const double G = 6.6743e-11;
const double M = 1.9891e30;

double t_start = 0.0;
double t_end = 10*365*24*3600.0;
double h = 100.0;

vector<double> m = {-4.0e11,5.0e11,1.0e4,-1.0e4}; //x,y,vx,vy

double a_x(double x, double y){
    return -G*M*x/pow(sqrt(pow(x,2)+pow(y,2)),3);
}
double a_y(double x, double y){
    return -G*M*y/pow(sqrt(pow(x,2)+pow(y,2)),3);
}

vector<double> rk4(double x0, double y0, double vx0, double vy0){
    double k1_x = h*vx0;
    double k1_y = h*vy0;
    double l1_x = h*a_x(x0, y0);
    double l1_y = h*a_y(x0, y0);
    double k2_x = h*(vx0 + l1_x/2.0);
    double k2_y = h*(vy0 + l1_y/2.0);
    double l2_x = h*a_x(x0 + k1_x/2.0, y0 + k1_y/2.0);
    double l2_y = h*a_y(x0 + k1_x/2.0, y0 + k1_y/2.0);
    double k3_x = h*(vx0 + l2_x/2.0);
    double k3_y = h*(vy0 + l2_y/2.0);
    double l3_x = h*a_x(x0 + k2_x/2.0, y0 + k2_y/2.0);
    double l3_y = h*a_y(x0 + k2_x/2.0, y0 + k2_y/2.0);
    double k4_x = h*(vx0 + l3_x);
    double k4_y = h*(vy0 + l3_y);
    double l4_x = h*a_x(x0 + k3_x, y0 + k3_y);
    double l4_y = h*a_y(x0 + k3_x, y0 + k3_y);
    double x = x0 + (k1_x + 2*k2_x + 2*k3_x + k4_x)/6.0;
    double vx = vx0 + (l1_x + 2*l2_x + 2*l3_x + l4_x)/6.0;
    double y = y0 + (k1_y + 2*k2_y + 2*k3_y + k4_y)/6.0;
    double vy = vy0 + (l1_y + 2*l2_y + 2*l3_y + l4_y)/6.0;
    vector<double> m_new = {x,y,vx,vy};
    return m_new;
}
int main()
{
    ofstream file;
    file.open("output.csv");
    for(double t=t_start;t<=t_end;t+=h){
        m = rk4(m[0],m[1],m[2],m[3]);
        file<<m[0]<<","<<m[1]<<","<<m[2]<<","<<m[3]<<'\n';
    }
    file.close();
    cout<<"Completed"<<endl;
    return 0;
}