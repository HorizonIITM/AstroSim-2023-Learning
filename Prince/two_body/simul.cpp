#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

const float G = 1;
const float step = 0.01;
int counter = 0;

/*
float Q_rsqrt( float number )
{
 long i;
 float x2, y;
 const float threehalfs = 1.5F;

 x2 = number * 0.5F;
 y  = number;
 i  = * ( long * ) &y;                       // to convert long without destroying the binary bits
 i  = 0x5f3759df - ( i >> 1 );               // moving the binary to right  - to reduce decimal by half
 y  = * ( float * ) &i;
 y  = y * ( threehalfs - ( x2 * y * y ) );   // Newton's iteration
  

 return y;
}
*/

float Q_rsqrt(float number){
    union {
    float    f;
    uint32_t i;
    } conv = { .f = number };
    conv.i  = 0x5f3759df - (conv.i >> 1);
    conv.f *= 1.5F - (number * 0.5F * conv.f * conv.f);
    return conv.f;}

double mag_sqr(const vector<float>& vector) {
    double magnitude = 0.0;
    
    for (size_t i = 0; i < vector.size(); i++) {
        magnitude += vector[i] * vector[i];
    }
    
    return magnitude;
};

class solver {
public:
    template <class func>
    void rk_4(func f, float t, vector<float>& x) {
        vector<float> k1(3), k2(3), k3(3), k4(3);
        
        for(int i=0; i<2; i++){
            k1[i] = step * f(t, x[i]);
            k2[i] = step * f(t + step/2.0, x[i] + step * k1[i]/2.0);
            k3[i] = step * f(t + step/2.0, x[i] + step * k2[i]/2.0);
            k4[i] = step * f(t + step, x[i] + step * k3[i]);
        
            x[i] += step * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6.0;

            counter++;
        }
        
    }
};

class gravitational_system {
public:
    float M, m, t;
    vector<float> q;
    vector<float> p;

    gravitational_system(float Mi, float mi, const vector<float> qi, const vector<float> pi, float ti)
        : M(Mi), m(mi), q(qi), p(pi), t(ti) {}

    void write_xy(ostream& outstream) {
        outstream << q[0] << "," << q[1] << endl;
    }
};

class q_deriv {
private:
    gravitational_system& sys;

public:
    q_deriv(gravitational_system& s) : sys(s) {}

    float operator()(float t, float q) { 
        if (counter == 0) {            
            return sys.p[0] / sys.m;        
        }
        else {
            return sys.p[1] / sys.m;
        }
    }
};

class p_deriv {
private:
    gravitational_system& sys;

public:
    p_deriv(gravitational_system& s) : sys(s) {}

    float operator()(float t, float p) {  
        if (counter == 0){
            return  -(G * sys.M * sys.m * sys.q[0]) * (pow(Q_rsqrt(mag_sqr(sys.q)), 3.0));
        }
        else{
            return  -(G * sys.M * sys.m * sys.q[1]) * (pow(Q_rsqrt(mag_sqr(sys.q)), 3.0));
        }

    }
};

class gravi_solver : private solver, public gravitational_system {

public:
    gravi_solver(float Mi, float mi,vector<float> qi, vector<float> pi, float ti):  gravitational_system(Mi,mi,qi,pi,ti){}

private:
    void rk_4() {
        counter = 0;
        solver::rk_4(q_deriv(*this), t, q);
        
        counter = 0;
        solver::rk_4(p_deriv(*this), t, p);
        
        t += step;
        
    }
public:
    void solve(float total_time, std::string filename = ""){
            bool write_flag = filename=="" ? false : true;
            std::ofstream my_file(filename);                       

            if(write_flag) my_file<<"x,y"<<std::endl;
            for(int i=0;i<total_time/step;i++){
                if(write_flag) write_xy(my_file);
                rk_4();
            }
            if(write_flag) write_xy(my_file);
            my_file.close();
        }
};

int main(int argc, char* argv[]) {

    float M = 100;
    float m = 1;
    float xi = 20, yi = 0, zi = 0, pxi = 0, pyi = sqrt(1.5 * G * M / xi), pzi = 0;
    float total_time = 20000;
    
    vector<float> q;
	q.push_back(xi);
	q.push_back(yi);

    vector<float> p;
	p.push_back(pxi);
	p.push_back(pyi);
	
    gravi_solver my_solver(M, m, q, p, 0);
    my_solver.solve(total_time, "draw.txt");
    
    return 0;
}
