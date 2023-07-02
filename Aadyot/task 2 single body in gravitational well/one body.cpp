#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;
const float G = 1;

float Q_rsqrt(float number){
    union {
    float    f;
    uint32_t i;
    } conv = { .f = number };
    conv.i  = 0x5f3759df - (conv.i >> 1);
    conv.f *= 1.5F - (number * 0.5F * conv.f * conv.f);
    return conv.f;
}

template <class func>
class solver{
    //a generic rk4 solver. takes function, current value, step
    //give [dx/dt = f(x,t), t, x at t, step] it finds x at (t+step)
    public:
        float x;
        float t;
        float step;
        func f;

        solver(){}  //bad coding here. some places wanted default constuctor. same problem in gravitational system
        solver(float xi, float t, float step, func f): x(xi), step(step),t(t), f(f) {}

        float next_step_rk4(){
            
            float k1 = step * f(t,x);  //other variables are constant!!
            float k2 = step * f(t + step/2.0 , x + step*k1/2.0);
            float k3 = step * f(t + step/2.0 , x + step*k2/2.0);
            float k4 = step * f(t + step , x + step*k3);

            x += (step/6.0)* (k1+ 2*k2 + 2*k3 + k4);
            t +=step;
            return x;
        }
};


class gravitational_system{
    //contains the info about the gravitational system
    public:
        float M,m,x,y,px,py,t;

        gravitational_system(){}   //bad coding here. some places wanted default constuctor. same problem in solver

        gravitational_system(float Mi, float mi,float xi, float yi, float pxi, float pyi,float ti)
        :M(Mi), m(mi), x(xi), y(yi), px(pxi), py(pyi),t(ti){}

        void print_xy(){
            cout<<x<<" "<<y<<endl;
        }

        void write_xy(ofstream& filestream){
            filestream<<x<<","<<y<<endl;
        }

        gravitational_system update_pos_mom(float xi, float yi, float pxi, float pyi, float ti){
            x = xi;
            y = yi;
            px = pxi;
            py = pyi;
            t = ti;
            return *this;
        }
};

namespace differential_equations{
    //differential equations governing the system
    class x_deriv : gravitational_system{
        public:
            x_deriv(){};
            x_deriv(gravitational_system s): gravitational_system(s){};
            float operator() (float t, float x){
                return px / m;
            }
    };
    class y_deriv : gravitational_system{
        public:
            y_deriv(){}
            y_deriv(gravitational_system s): gravitational_system(s){}
            float operator() (float t, float y){
                return py / m;
            }
    };
    class px_deriv : gravitational_system{
        public:
            px_deriv(){}
            px_deriv(gravitational_system s): gravitational_system(s){}
            float operator() (float t, float px){
                return - (G * M * m * x) *(pow (Q_rsqrt(x*x + y*y),3.0));
            }
    };
    class py_deriv : gravitational_system{
        public:
            py_deriv(){}
            py_deriv(gravitational_system s): gravitational_system(s){}
            float operator() (float t, float py){
                return - (G * M * m * y) * (pow (Q_rsqrt(x*x + y*y),3.0));
            }
    };
}




class gravitational_solver{
    public:
        //in future codes, this can be made more general by taking any system as template arg
        float step;

        gravitational_system current_system;
        solver<differential_equations::x_deriv> x_solver;
        solver<differential_equations::y_deriv> y_solver;
        solver<differential_equations::px_deriv> px_solver;
        solver<differential_equations::py_deriv> py_solver;
        
        gravitational_solver(gravitational_system s, float step){
            this->step = step;
            current_system = s;
            // x_solver = solver(s.x, s.t, step, differential_equations::x_deriv(s));
            // y_solver = solver(s.y, s.t, step, differential_equations::y_deriv(s));
            // px_solver = solver(s.px, s.t, step, differential_equations::px_deriv(s));
            // py_solver = solver(s.py, s.t, step, differential_equations::py_deriv(s));
        }

        gravitational_system next_step_all_rk4(){

            // i want to avoid initialization here and i want to do it just once in the constructor
            x_solver = solver(current_system.x, current_system.t, step, differential_equations::x_deriv(current_system));
            y_solver = solver(current_system.y, current_system.t, step, differential_equations::y_deriv(current_system));
            px_solver = solver(current_system.px, current_system.t, step, differential_equations::px_deriv(current_system));
            py_solver = solver(current_system.py, current_system.t, step, differential_equations::py_deriv(current_system));
            
            current_system= current_system.update_pos_mom(
                x_solver.next_step_rk4(),
                y_solver.next_step_rk4(),
                px_solver.next_step_rk4(),
                py_solver.next_step_rk4(),
                current_system.t+step
            );
            return current_system;
        }
};


int main(int argc, char *argv[]) {
    //can i take all of these as command line args?
    float central_mass = 100;
    float moving_mass = 1;
    float xi = 20, yi = 0, pxi = 0, pyi = sqrt(G*central_mass/xi);
    float step = 0.01;
    float total_time = 9000;

    if(argc==2) total_time = stoi(argv[1]);
    cout<<total_time<<endl;


    gravitational_system my_system(central_mass, moving_mass, xi, yi, pxi, pyi, 0);
    
    gravitational_solver my_solver(my_system, step);

    ofstream my_file("try1.txt");

    my_file<<"x,y"<<endl;
    for(int i=0;i<total_time/step;i++){
         my_system.write_xy(my_file);
         my_system = my_solver.next_step_all_rk4();
    }
    my_system.write_xy(my_file);
    
    my_file.close();

    return 0;
}