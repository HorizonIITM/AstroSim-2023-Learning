#include <iostream>
#include <cmath>
#include <fstream>

#define ll long double

using namespace std;

const ll G = 6.6743 * pow(10,-11);

const ll M = pow(10, 14);                          //mass in kg
pair <ll, ll> v = make_pair(0, 10);   //initial velocity (vx, vy)
pair <ll, ll> r = make_pair(66.74, 0);     //initial position (x, y)

const ll h = 0.001;                         //step size



ll ax(ll x, ll y){        
    return -G * M * x/pow(sqrt(x*x + y*y), 3);
}

ll ay(ll x, ll y){        
    return -G * M * y/pow(sqrt(x*x + y*y), 3);
}

void RungeX(ll *x, ll *y, ll *v){
    ll k1 = h * *v;
    ll l1 = h * ax(*x, *y);
    ll k2 = h * (*v + l1/2);
    ll l2 = h * ax(*x + k1/2, *y);
    ll k3 = h * (*v + l2/2);
    ll l3 = h * ax(*x + k2/2, *y);
    ll k4 = h * (*v + l3);
    ll l4 = h * ax(*x + k3, *y);
    *v += (l1 + 2*l2 + 2*l3 + l4)/6;
    *x += (k1 + 2*k2 + 2*k3 + k4)/6;
}

void RungeY(ll *x, ll *y, ll *v){
    ll k1 = h * *v;
    ll l1 = h * ay(*x, *y);
    ll k2 = h * (*v + l1/2);
    ll l2 = h * ay(*x, *y + k1/2);
    ll k3 = h * (*v + l2/2);
    ll l3 = h * ay(*x, *y + k2/2);
    ll k4 = h * (*v + l3);
    ll l4 = h * ay(*x , *y + k3);
    *v += (l1 + 2*l2 + 2*l3 + l4)/6;
    *y += (k1 + 2*k2 + 2*k3 + k4)/6;
}

int main(){
    ofstream file("2BodySystem.txt");
    file << r.first << ' ' << r.second << endl;
    
    for(int i = 0; i < 50000; i++){                            // number of points
        RungeX(&r.first, &r.second, &v.first);
        RungeY(&r.first, &r.second, &v.second);
               
        file << r.first << ' ' << r.second << endl;  
    }

    file.close();
}
