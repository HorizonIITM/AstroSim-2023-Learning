#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
const float G = 1;
const float step = 0.01;

float Q_rsqrt(float number){
    union {
    float    f;
    uint32_t i;
    } conv = { .f = number };
    conv.i  = 0x5f3759df - (conv.i >> 1);
    conv.f *= 1.5F - (number * 0.5F * conv.f * conv.f);
    return conv.f;
}

class solver{
    //a generic rk4 solver. takes function, current value, step
    //give [dx/dt = f(x,t), t, x at t, step] it finds x at (t+step)
    public:
        template <class func>
        void next_step_rk4_set(float& x, float t, func f){
            float k1 = step * f(t,x);  //other variables are constant!!
            float k2 = step * f(t + step/2.0 , x + step*k1/2.0);
            float k3 = step * f(t + step/2.0 , x + step*k2/2.0);
            float k4 = step * f(t + step , x + step*k3);

            x += (step/6.0)* (k1+ 2*k2 + 2*k3 + k4);
            t +=step;
        }
};

class gravitational_system{
    //contains the info about the gravitational system
    public:
        float M,m,x,y,px,py,t;

        gravitational_system(float Mi, float mi,float xi, float yi, float pxi, float pyi,float ti)
        :M(Mi), m(mi), x(xi), y(yi), px(pxi), py(pyi),t(ti){}

        void write_xy(std::ostream& outstream){
            outstream<<x<<","<<y<<std::endl;
        }
};

namespace differential_equations{
    //differential equations governing the system
    class x_deriv{
        private:
            gravitational_system& sys;
        public:
            x_deriv(gravitational_system s): sys(s){}
            float operator() (float t, float x){
                return sys.px / sys.m;
            }
    };
    class y_deriv{
        private:
            gravitational_system& sys;
        public:
            y_deriv(gravitational_system s): sys(s){}
            float operator() (float t, float y){
                return sys.py / sys.m;
            }
    };
    class px_deriv{
        private:
            gravitational_system& sys;
        public:
            px_deriv(gravitational_system s): sys(s){}
            float operator() (float t, float px){
                return - (G * sys.M * sys.m * sys.x) *(pow (Q_rsqrt(sys.x*sys.x + sys.y*sys.y),3.0));
            }
    };
    class py_deriv{
        private:
            gravitational_system& sys;
        public:
            py_deriv(gravitational_system s): sys(s){}
            float operator() (float t, float py){
                return - (G * sys.M * sys.m * sys.y) * (pow (Q_rsqrt(sys.x*sys.x + sys.y*sys.y),3.0));
            }
    };
}

class gravitational_solver : private solver, public gravitational_system{
    public:
        gravitational_solver(float Mi, float mi,float xi, float yi, float pxi, float pyi,float ti):  gravitational_system(Mi,mi,xi,yi,pxi,pyi,ti){}
    
    private:
        void next_step_all_rk4_set(){ 
            next_step_rk4_set(x, t, differential_equations::x_deriv(*this));
            next_step_rk4_set(y, t, differential_equations::y_deriv(*this));
            next_step_rk4_set(px, t, differential_equations::px_deriv(*this));
            next_step_rk4_set(py, t, differential_equations::py_deriv(*this));                         
            t+=step;
        }
    public: 
        void solve(float total_time, std::string filename = ""){
            bool write_flag = filename=="" ? false : true;
            std::ofstream my_file(filename);                       

            if(write_flag) my_file<<"x,y"<<std::endl;
            for(int i=0;i<total_time/step;i++){
                if(write_flag) write_xy(my_file);
                next_step_all_rk4_set();
            }
            if(write_flag) write_xy(my_file);
            my_file.close();
        }
};


int main(int argc, char *argv[]) {
    //can i take all of these as command line args?
    float central_mass = 100;
    float moving_mass = 1;
    float xi = 20, yi = 0, pxi = 0, pyi = sqrt(1.5*G*central_mass/xi);
    float total_time = 20000;
    gravitational_solver my_solver(central_mass, moving_mass, xi, yi, pxi, pyi, 0);
    my_solver.solve(total_time, "try3.txt");
    return 0;
}