#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;
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
        /* template <class func>
        float next_step_rk4(float x, float t, func f){
            
            float k1 = step * f(t,x);  //other variables are constant!!
            float k2 = step * f(t + step/2.0 , x + step*k1/2.0);
            float k3 = step * f(t + step/2.0 , x + step*k2/2.0);
            float k4 = step * f(t + step , x + step*k3);

            x += (step/6.0)* (k1+ 2*k2 + 2*k3 + k4);
            t +=step;
            return x;
        }*/
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

        void write_xy(ofstream& outstream){
            outstream<<x<<","<<y<<endl;
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




class gravitational_solver : private solver{
    private:
        gravitational_system& current_system;
        // differential_equations::x_deriv xp{current_system};
        // differential_equations::y_deriv yp{current_system};
        // differential_equations::px_deriv pxp{current_system};
        // differential_equations::py_deriv pyp{current_system};

    public:
        //in future codes, this can be made more general by taking any system as template arg
        gravitational_solver(gravitational_system& s): current_system(s){}
    
    private:
        void next_step_all_rk4_set(){  
            next_step_rk4_set(current_system.x, current_system.t, differential_equations::x_deriv{current_system});
            next_step_rk4_set(current_system.y, current_system.t, differential_equations::y_deriv{current_system});
            next_step_rk4_set(current_system.px, current_system.t, differential_equations::px_deriv{current_system});
            next_step_rk4_set(current_system.py, current_system.t, differential_equations::py_deriv{current_system});                         
            current_system.t+=step;
            // next_step_rk4_set(current_system.x, current_system.t, xp);
            // next_step_rk4_set(current_system.y, current_system.t, yp);
            // next_step_rk4_set(current_system.px, current_system.t, pxp);
            // next_step_rk4_set(current_system.py, current_system.t, pyp);                         
            // current_system.t+=step;
        }
    public: 
        void solve(float total_time, string filename = ""){
            bool write_flag = true;
            ofstream my_file(filename);
            if(filename=="") {
                write_flag = false;
                my_file.close();
            }
                       

            if(write_flag) my_file<<"x,y"<<endl;
            for(int i=0;i<total_time/step;i++){
                if(write_flag) current_system.write_xy(my_file);
                next_step_all_rk4_set();
            }
            if(write_flag) current_system.write_xy(my_file);
            if(write_flag) my_file.close();
        }


};


int main(int argc, char *argv[]) {
    //can i take all of these as command line args?
    float central_mass = 100;
    float moving_mass = 1;
    float xi = 20, yi = 0, pxi = 0, pyi = sqrt(1.5*G*central_mass/xi);
    float total_time = 20000;

    gravitational_system my_system(central_mass, moving_mass, xi, yi, pxi, pyi, 0);
    
    gravitational_solver my_solver(my_system);

    my_solver.solve(total_time, "try3.txt");

    return 0;
}