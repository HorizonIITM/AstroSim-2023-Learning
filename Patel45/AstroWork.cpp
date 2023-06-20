#include<bits/stdc++.h>
using namespace std;

//Universal Constants
// const double G=6.67430e-11;
const double G=1;

//Useful function
double invSqrt(double);

//Using Cartesian Coordinates
class celestialBody;
class movingCelestialBody;

class celestialBody{
    //Star/Sun in this case is assumed to be stationary due to its large size.
    protected: 
        string name;
        double mass;
        vector<double> position;
    public:
        //Constructors 
        celestialBody(string s,double m,vector<double> p){
            name=s;
            mass=m;
            position=p;
        }
        //Friend Function
        friend vector<double> Force(movingCelestialBody,celestialBody);
        friend void move(movingCelestialBody&,celestialBody,double);
        friend void basicDEsolver(vector<double>,vector<double>&,double);
        //Show State
        void currentState(){
            cout<<name<<":\n";
            cout<<setprecision(6)<<"Mass : "<<mass<<"\n";
            cout<<setprecision(6)<<"Position : ("<<position[0]<<","<<position[1]<<","<<position[2]<<")\n";
            cout<<"\n";
        }
};

class movingCelestialBody : public celestialBody{
    //Bodies in Motion like Planets, Moon, Asteriods, etc.
    protected:
        vector<double> velocity; 
    public:
        //Constructor 
        movingCelestialBody(string s,double m, vector<double> p, vector<double> v) : celestialBody(s,m,p){
            velocity=v;
        }
        //Friend Function
        friend vector<double> Force(movingCelestialBody,celestialBody);
        friend void move(movingCelestialBody&,celestialBody,double);
        friend void basicDEsolver(vector<double>,vector<double>&,double);
        //Show State
        void currentState(){
            cout<<name<<":\n";
            cout<<setprecision(6)<<"Mass : "<<mass<<"\n";
            cout<<setprecision(6)<<"Position : ("<<position[0]<<","<<position[1]<<","<<position[2]<<")\n";
            cout<<setprecision(6)<<"Velocity : ("<<velocity[0]<<","<<velocity[1]<<","<<velocity[2]<<")\n";
            cout<<"\n";
        }
        void outToFile(ofstream& outdata){
            outdata<<setprecision(10)<<position[0]<<" "<<position[1]<<" "<<position[2]<<"\n";
        }
};

vector<double> Force(movingCelestialBody b1, celestialBody b2){
    vector<double> force(3);
    double distSq=pow((b1.position[0]-b2.position[0]),2)+pow((b1.position[1]-b2.position[1]),2)+pow((b1.position[2]-b2.position[2]),2);
    for(int i=0; i<3; i++){
        force[i]=G*b1.mass*b2.mass*invSqrt(distSq)*(b2.position[i]-b1.position[i])/distSq;
    }
    return force;
}

//Use the desired function.
double invSqrt(double x){
    return 1.0/sqrt(x);
}

void basicDEsolver(vector<double> derivative, vector<double> &current_value, double step_size){
    for(int d=0; d<3; d++){
        current_value[d]+=derivative[d]*step_size;
    }
}
void move(movingCelestialBody &b1, celestialBody b2, double step_size){
    // Get all the forces acting on the body
    //Here we have only 1 force 
    vector<double> acceleration;
    for(auto f:Force(b1,b2)) acceleration.push_back(f/b1.mass);
    basicDEsolver(b1.velocity,b1.position,step_size);
    basicDEsolver(acceleration,b1.velocity,step_size);
}

//For Output to Terminal
void output(vector<double> v){
    cout<<"("<<v[0]<<","<<v[1]<<","<<v[2]<<")\n";
}

//For output to File
int getOpen(ofstream& fileOut)
{
   string name;
   cout << "\nEnter a file name: ";
   getline(cin,name);
   fileOut.open(name.c_str());   // open the file
   if (fileOut.fail())   // check for successful open
   {
      cout << "Cannot open the file" << endl;
      exit(1);
   }
   else
      return 1;
}

int main(){
    ofstream outdata;
    getOpen(outdata);
    
    //Actual
    //SI Units Used
    // celestialBody Sun("Sun",1.9891e30,{0,0,0});
    // movingCelestialBody Earth("Earth",5.972e24,{1.5e11,0,0},{0,,0});
    
    //Star
    celestialBody Sun("Sun",10,{0,0,0});
    
    //UCM - Earth
    // movingCelestialBody Earth("Earth",1,{10,0,0},{0,1,0});
    
    // Parabolic Path - Earth
    movingCelestialBody Earth("Earth",1,{10,0,0},{0,1.41421356237,0});
    
    cout<<"Initial State \n";
    Sun.currentState();
    Earth.currentState();
    
    cout<<"Simulation Begins\n\n";
    for(int i=0; i<=30000; i++){
        Earth.outToFile(outdata);
        if(i%1000==0){
            cout<<"Step "<<i<<":\n";
            cout<<"Force Acting ";
            output(Force(Earth,Sun));
            Earth.currentState();   
            cout<<"\n"; 
        }
        move(Earth,Sun,1e-2);
    }
    outdata.close();
    cout<<"Simulation Ends\n";
    return 0;
}