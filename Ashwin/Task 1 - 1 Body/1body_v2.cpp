#include <bits/stdc++.h>
#include <vector>
using namespace std;
const double G = 6.6743e-11;
const double M = 1.9891e30;

class Body {
    public:
        double mass;
        vector<double> position;
        vector<double> momentum;
        Body(double m,double x,double y,double px,double py){
            mass = m;
            position = {x,y};
            momentum = {px,py};
        }
};

class solver { // this solver is meant for the one body case
    private:
        double px_deriv(double m, double x, double y){
            return -(G*M*m*x)/pow(x*x+y*y,1.5);
        }
        double py_deriv(double m, double x,double y){
            return -(G*M*m*y)/pow(x*x+y*y,1.5);
        }
        double x_deriv(double m, double px){
            return px/m;
        }
        double y_deriv(double m, double py){
            return py/m;
        }
        vector<double> k1(double h, Body body){
            double m = body.mass;
            double x1 = body.position[0];
            double y1 = body.position[1];
            double px1 = body.momentum[0];
            double py1 = body.momentum[1];
            vector<double> k1_all = {h*x_deriv(m,px1),h*y_deriv(m,py1),h*px_deriv(m,x1,y1),h*py_deriv(m,x1,y1)}; //{k1x,k1y,k1px,k1py}
            return k1_all;
        }
        vector<double> k2(double h, Body body, vector<double> k1_all){
            double m = body.mass;
            double x2 = body.position[0] + k1_all[0]/2;
            double y2 = body.position[1] + k1_all[1]/2;
            double px2 = body.momentum[0] + k1_all[2]/2;
            double py2 = body.momentum[1] + k1_all[3]/2;
            vector<double> k2_all = {h*x_deriv(m,px2),h*y_deriv(m,py2),h*px_deriv(m,x2,y2),h*py_deriv(m,x2,y2)}; //{k2x,k2y,k2px,k2py}
            return k2_all;
        }
        vector<double> k3(double h, Body body, vector<double> k2_all){
            double m = body.mass;
            double x3 = body.position[0] + k2_all[0]/2;
            double y3 = body.position[1] + k2_all[1]/2;
            double px3 = body.momentum[0] + k2_all[2]/2;
            double py3 = body.momentum[1] + k2_all[3]/2;
            vector<double> k3_all = {h*x_deriv(m,px3),h*y_deriv(m,py3),h*px_deriv(m,x3,y3),h*py_deriv(m,x3,y3)}; //{k3x,k3y,k3px,k3py}
            return k3_all;
        }
        vector<double> k4(double h, Body body, vector<double> k3_all){
            double m = body.mass;
            double x3 = body.position[0] + k3_all[0];
            double y3 = body.position[1] + k3_all[1];
            double px3 = body.momentum[0] + k3_all[2];
            double py3 = body.momentum[1] + k3_all[3];
            vector<double> k4_all = {h*x_deriv(m,px3),h*y_deriv(m,py3),h*px_deriv(m,x3,y3),h*py_deriv(m,x3,y3)}; //{k4x,k4y,k4px,k4py}
            return k4_all;
        }
    public:
        Body rk4(Body body, double h){
            vector<double> k1_body = k1(h,body);
            vector<double> k2_body = k2(h,body,k1_body);
            vector<double> k3_body = k2(h,body,k2_body);
            vector<double> k4_body = k2(h,body,k3_body);
            body.position[0] += (k1_body[0] + 2*k2_body[0] + 2*k3_body[0] + k4_body[0])/6;
            body.position[1] += (k1_body[1] + 2*k2_body[1] + 2*k3_body[1] + k4_body[1])/6;
            body.momentum[0] += (k1_body[2] + 2*k2_body[2] + 2*k3_body[2] + k4_body[2])/6;
            body.momentum[1] += (k1_body[3] + 2*k2_body[3] + 2*k3_body[3] + k4_body[3])/6;
            return body;
        }
};

int main(){
    Body earth = Body(5.972e24,2.49e11,0,0,1.2213e4*5.972e24);
    double t_start = 0;
    double t_end = 365*24*3600;
    double h = 6;
    solver gravity_solver;
    ofstream file;
    file.open("output_ell.csv");
    for(double t=t_start;t<=t_end;t+=h){
        earth = gravity_solver.rk4(earth,h);
        file<<earth.position[0]<<","<<earth.position[1]<<","<<earth.momentum[0]<<","<<earth.momentum[1]<<'\n';
    }
    file.close();
    cout<<"Completed"<<endl;
    return 0;
}
