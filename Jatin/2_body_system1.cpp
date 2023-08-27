#include <iostream>
#include <cmath>
#include <fstream>

#define ll long double

using namespace std;

const ll G = 6.6743 * pow(10, -11);

const ll h = 0.001;

class Gsystem
{
public:
    ll M = pow(10, 14), m = 1, t = 0;     // mass in kg
    pair<ll, ll> p = make_pair(0, 10);    // initial momentum (mx, my)
    pair<ll, ll> r = make_pair(66.74, 0); // initial position (x, y)
} S;

class Vx
{
private:
    Gsystem &sys;

public:
    Vx(Gsystem a) : sys(a){};
    ll operator()(ll x, ll t)
    {
        return S.p.first / S.m;
    }
};

class Vy
{
private:
    Gsystem &sys;

public:
    Vy(Gsystem a) : sys(a){};
    ll operator()(ll x, ll t)
    {
        return S.p.second / S.m;
    }
};

class Fx
{
private:
    Gsystem &sys;

public:
    Fx(Gsystem a) : sys(a){};
    ll operator()(ll x, ll t)
    {
        return -G * S.M * S.m * S.r.first / pow(sqrt(S.r.first * S.r.first + S.r.second * S.r.second), 3);
    }
};

class Fy
{
private:
    Gsystem &sys;

public:
    Fy(Gsystem a) : sys(a){};
    ll operator()(ll x, ll t)
    {
        return -G * S.M * S.m * S.r.second / pow(sqrt(S.r.first * S.r.first + S.r.second * S.r.second), 3);
    }
};

class Differential_Solver
{
public:
    template <class func>
    void runge(ll &x, ll t, func f)
    {
        float k1 = f(x, t);
        float k2 = f(x + h * k1 / 2, t + h / 2);
        float k3 = f(x + h * k2 / 2, t + h / 2);
        float k4 = f(x + h * k3, t + h);
        t += h;

        x += (k1 + 2 * k2 + 2 * k3 + k4) * h / 6;
    }
} solve;

int main()
{
    ofstream file("2BodySystem.txt");
    file << S.r.first << ' ' << S.r.second << endl;

    for (int i = 0; i < 50000; i++)
    {
        solve.runge(S.r.first, S.t, Vx{S});
        solve.runge(S.r.second, S.t, Vy{S});
        solve.runge(S.p.first, S.t, Fx{S});
        solve.runge(S.p.second, S.t, Fy{S});

        file << S.r.first << ' ' << S.r.second << endl;
    }

    file.close();
}