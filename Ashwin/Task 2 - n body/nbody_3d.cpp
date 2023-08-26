#include <bits/stdc++.h>
#include <vector>
using namespace std;
const double G = 1;

class Body {
    public:
        double mass;
        vector<double> position;
        vector<double> momentum;
        Body(double m,double x,double y,double z,double px,double py,double pz){
            mass = m;
            position = {x,y,z};
            momentum = {px,py,pz};
        }
};

class solver {
    private:
        vector<double> px_deriv(vector<Body> o, int n){
            vector<double> pxd(n,0);
            for (int i=0;i<n-1;i++){
                for (int j=i+1;j<n;j++){
                    double xr = o[j].position[0] - o[i].position[0];
                    double yr = o[j].position[1] - o[i].position[1];
                    double zr = o[j].position[2] - o[i].position[2];
                    pxd[i]+=(G*o[i].mass*o[j].mass*xr)/pow(xr*xr+yr*yr+zr*zr,1.5);
                    pxd[j]-=(G*o[i].mass*o[j].mass*xr)/pow(xr*xr+yr*yr+zr*zr,1.5);
                }
            }
            return pxd;
        }
        vector<double> py_deriv(vector<Body> o, int n){
            vector<double> pyd(n,0);
            for (int i=0;i<n-1;i++){
                for (int j=i+1;j<n;j++){
                    double xr = o[j].position[0] - o[i].position[0];
                    double yr = o[j].position[1] - o[i].position[1];
                    double zr = o[j].position[2] - o[i].position[2];
                    pyd[i]+=(G*o[i].mass*o[j].mass*yr)/pow(xr*xr+yr*yr+zr*zr,1.5);
                    pyd[j]-=(G*o[i].mass*o[j].mass*yr)/pow(xr*xr+yr*yr+zr*zr,1.5);
                }
            }
            return pyd;
        }
        vector<double> pz_deriv(vector<Body> o, int n){
            vector<double> pzd(n,0);
            for (int i=0;i<n-1;i++){
                for (int j=i+1;j<n;j++){
                    double xr = o[j].position[0] - o[i].position[0];
                    double yr = o[j].position[1] - o[i].position[1];
                    double zr = o[j].position[2] - o[i].position[2];
                    pzd[i]+=(G*o[i].mass*o[j].mass*zr)/pow(xr*xr+yr*yr+zr*zr,1.5);
                    pzd[j]-=(G*o[i].mass*o[j].mass*zr)/pow(xr*xr+yr*yr+zr*zr,1.5);
                }
            }
            return pzd;
        }
        vector<double> x_deriv(vector<Body> o, int n){
            vector<double> xd(n,0);
            for (int i=0;i<n;i++){
                xd[i]=o[i].momentum[0]/o[i].mass;
            }
            return xd;
        }
        vector<double> y_deriv(vector<Body> o, int n){
            vector<double> yd(n,0);
            for (int i=0;i<n;i++){
                yd[i]=o[i].momentum[1]/o[i].mass;
            }
            return yd;
        }
        vector<double> z_deriv(vector<Body> o, int n){
            vector<double> zd(n,0);
            for (int i=0;i<n;i++){
                zd[i]=o[i].momentum[2]/o[i].mass;
            }
            return zd;
        }
        vector<Body> k1(vector<Body> o, int n, double h){
            vector<double> xd1 = x_deriv(o,n);
            vector<double> yd1 = y_deriv(o,n);
            vector<double> zd1 = z_deriv(o,n);
            vector<double> pxd1 = px_deriv(o,n);
            vector<double> pyd1 = py_deriv(o,n);
            vector<double> pzd1 = pz_deriv(o,n);

            vector<Body> k1_body = o;
            for (int i=0;i<n;i++){
                k1_body[i].position[0] += h*xd1[i]/2; //x2
                k1_body[i].position[1] += h*yd1[i]/2; //y2
                k1_body[i].position[2] += h*zd1[i]/2; //z2
                k1_body[i].momentum[0] += h*pxd1[i]/2; //px2
                k1_body[i].momentum[1] += h*pyd1[i]/2; //py2
                k1_body[i].momentum[2] += h*pzd1[i]/2; //py2
            }
            return k1_body;
        }
        vector<Body> k2(vector<Body> o, int n, double h, vector<Body> k1_body){
            vector<double> xd2 = x_deriv(k1_body,n);
            vector<double> yd2 = y_deriv(k1_body,n);
            vector<double> zd2 = z_deriv(k1_body,n);
            vector<double> pxd2 = px_deriv(k1_body,n);
            vector<double> pyd2 = py_deriv(k1_body,n);
            vector<double> pzd2 = pz_deriv(k1_body,n);
            vector<Body> k2_body = o;
            for (int i=0;i<n;i++){
                k2_body[i].position[0] += h*xd2[i]/2; //x3
                k2_body[i].position[1] += h*yd2[i]/2; //y3
                k2_body[i].position[2] += h*zd2[i]/2; //z3
                k2_body[i].momentum[0] += h*pxd2[i]/2; //px3
                k2_body[i].momentum[1] += h*pyd2[i]/2; //py3
                k2_body[i].momentum[2] += h*pzd2[i]/2; //pz3
            }
            return k2_body;
        }
        vector<Body> k3(vector<Body> o, int n, double h, vector<Body> k2_body){
            vector<double> xd3 = x_deriv(k2_body,n);
            vector<double> yd3 = y_deriv(k2_body,n);
            vector<double> zd3 = z_deriv(k2_body,n);
            vector<double> pxd3 = px_deriv(k2_body,n);
            vector<double> pyd3 = py_deriv(k2_body,n);
            vector<double> pzd3 = pz_deriv(k2_body,n);
            vector<Body> k3_body = o;
            for (int i=0;i<n;i++){
                k3_body[i].position[0] += h*xd3[i]; //x4
                k3_body[i].position[1] += h*yd3[i]; //y4
                k3_body[i].position[2] += h*zd3[i]; //z4
                k3_body[i].momentum[0] += h*pxd3[i]; //px4
                k3_body[i].momentum[1] += h*pyd3[i]; //py4
                k3_body[i].momentum[2] += h*pzd3[i]; //pz4
            }
            return k3_body;
        }
        vector<Body> k4(vector<Body> o, int n, double h, vector<Body> k3_body){
            vector<double> xd4 = x_deriv(k3_body,n);
            vector<double> yd4 = y_deriv(k3_body,n);
            vector<double> zd4 = z_deriv(k3_body,n);
            vector<double> pxd4 = px_deriv(k3_body,n);
            vector<double> pyd4 = py_deriv(k3_body,n);
            vector<double> pzd4 = pz_deriv(k3_body,n);
            vector<Body> k4_body = o;
            for (int i=0;i<n;i++){
                k4_body[i].position[0] += h*xd4[i]; 
                k4_body[i].position[1] += h*yd4[i]; 
                k4_body[i].position[2] += h*zd4[i];
                k4_body[i].momentum[0] += h*pxd4[i];
                k4_body[i].momentum[1] += h*pyd4[i];
                k4_body[i].momentum[2] += h*pzd4[i];
            }
            return k4_body;
        }
    public:
        vector<Body> rk4(vector<Body> o, int n, double h){
            vector<Body> k1_body = k1(o,n,h);
            vector<Body> k2_body = k2(o,n,h,k1_body);
            vector<Body> k3_body = k3(o,n,h,k1_body);
            vector<Body> k4_body = k4(o,n,h,k1_body);
            for (int i=0;i<n;i++){
                o[i].position[0] += (2*k1_body[i].position[0]+4*k2_body[i].position[0]+2*k3_body[i].position[0]+k4_body[i].position[0]-9*o[i].position[0])/6;
                o[i].position[1] += (2*k1_body[i].position[1]+4*k2_body[i].position[1]+2*k3_body[i].position[1]+k4_body[i].position[1]-9*o[i].position[1])/6;
                o[i].position[2] += (2*k1_body[i].position[2]+4*k2_body[i].position[2]+2*k3_body[i].position[2]+k4_body[i].position[2]-9*o[i].position[2])/6;
                o[i].momentum[0] += (2*k1_body[i].momentum[0]+4*k2_body[i].momentum[0]+2*k3_body[i].momentum[0]+k4_body[i].momentum[0]-9*o[i].momentum[0])/6;
                o[i].momentum[1] += (2*k1_body[i].momentum[1]+4*k2_body[i].momentum[1]+2*k3_body[i].momentum[1]+k4_body[i].momentum[1]-9*o[i].momentum[1])/6;
                o[i].momentum[2] += (2*k1_body[i].momentum[2]+4*k2_body[i].momentum[2]+2*k3_body[i].momentum[2]+k4_body[i].momentum[2]-9*o[i].momentum[2])/6; 
            }
            return o;
        }
};

int main(){
    vector<Body> o ={Body(1,1,0,0,0,0,1),Body(1,0,1,0,0,0,-1)};
    double t_start = 0;
    double t_end = 10;
    double h = 0.001;
    solver gravity_solver;
    ofstream file1;
    file1.open("foo1.csv");
    ofstream file2;
    file2.open("foo2.csv");
    file1<<"x,y,z"<<'\n';
    file2<<"x,y,z"<<'\n';
    for(double t=t_start;t<=t_end;t+=h){
        o = gravity_solver.rk4(o,2,h);
        file1<<o[0].position[0]<<","<<o[0].position[1]<<","<<o[0].position[2]<<'\n';
        file2<<o[1].position[0]<<","<<o[1].position[1]<<","<<o[1].position[2]<<'\n';
    }   
    file1.close();
    file2.close();
    cout<<"Completed"<<endl;
    return 0;
}