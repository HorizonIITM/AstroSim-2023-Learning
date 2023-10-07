#include <iostream>
#include <cmath>
#include <fstream>

#define ld long double

using namespace std;

const ld G = 6.6743 * pow(10, -11);     // gravitational constant

const ld h = 0.001;      //step size

const int n = 2;          //number of bodies

class Gsystem
{
public:
    ld M = 1, t = 0;                            //mass and time
    pair<ld, ld> p = make_pair(0, 0);           //momentum
    pair<ld, ld> r = make_pair(0, 0);           //position
}S[n];



class Vx
{
private:
    Gsystem &sys;

public:
    Vx(Gsystem a) : sys(a){};
    ld operator()(ld x, ld t, int k)
    {
        return S[k].p.first / S[k].M;
    }
};

class Vy
{
private:
    Gsystem &sys;

public:
    Vy(Gsystem a) : sys(a){};
    ld operator()(ld x, ld t, int k)
    {
        return S[k].p.second / S[k].M;
    }
};

class Fx
{
private:
    Gsystem &sys;

public:
    Fx(Gsystem a) : sys(a){};
    ld Fnet = 0, x = 0, y = 0, z = 0;
    ld operator()(ld x, ld t, int k)               //k is object index
    {
        for(int i = 0; i < n; i++){
            if(i == k) continue;
            x = S[k].r.first - S[i].r.first;
            y = S[k].r.second - S[i].r.second;

            Fnet += -G * S[i].M * S[k].M * x / pow(sqrt(x * x + y * y), 3);

        }

        return Fnet;
    }
};

class Fy
{
private:
    Gsystem &sys;

public:
    Fy(Gsystem a) : sys(a){};
    ld Fnet = 0, x = 0, y = 0, z = 0;
    ld operator()(ld x, ld t, int k)               //k is object index
    {
        for(int i = 0; i < n; i++){
            if(i == k) continue;
            x = S[k].r.first - S[i].r.first;
            y = S[k].r.second - S[i].r.second;

            Fnet += -G * S[i].M * S[k].M * y / pow(sqrt(x * x + y * y), 3);

        }

        return Fnet;
    }
};

class Differential_Solver
{
public:
    template <class func>
    void runge(ld &x, ld t, func f, int k)
    {
        float k1 = f(x, t, k);
        float k2 = f(x + h * k1 / 2, t + h / 2, k);
        float k3 = f(x + h * k2 / 2, t + h / 2, k);
        float k4 = f(x + h * k3, t + h, k);
        t += h;

        x += (k1 + 2 * k2 + 2 * k3 + k4) * h / 6;
    }
} solve;



int main()
{

    S[0].r = {-13, 0};
    S[0].p = {0, -1000};
    S[0].M = 1000000000;



    S[1].r = {13, 0};
    S[1].p = {0, 1000};
    S[1].M = 1000000000;





    ofstream file("NBodySystem.txt");

    for(int i = 0; i < n; i++){
        file << S[i].r.first << ',' << S[i].r.second ;
        if(i != n - 1) file << ",";
    }

    file << endl;

    for (int i = 0; i < 1000000; i++)
    {
        for(int k = 0; k < n; k++){
            solve.runge(S[k].r.first, S[k].t, Vx{S[k]}, k);
            solve.runge(S[k].r.second, S[k].t, Vy{S[k]}, k);
            solve.runge(S[k].p.first, S[k].t, Fx{S[k]}, k);
            solve.runge(S[k].p.second, S[k].t, Fy{S[k]}, k);
        }

        // solve.runge(S.r.first, S.t, Vx{S});
        // solve.runge(S.r.second, S.t, Vy{S});
        // solve.runge(S.p.first, S.t, Fx{S});
        // solve.runge(S.p.second, S.t, Fy{S});


        for(int k = 0; k < n; k++){
            file << S[k].r.first << ',' << S[k].r.second;
            if(k != n - 1) file << ",";
        }
        file << endl;
        
    }

    file.close();
}



