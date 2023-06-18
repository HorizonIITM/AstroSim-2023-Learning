template <class d>
class solver{
    public:
        solver(){}   //how to get rid of this
        solver(int a, int b, d f){}
};

class egg{
    public:
        float operator()(){
            return 1;
        }
};

class gsol{
    solver<egg> x_solver;
    public:
        gsol(){
            x_solver = solver<egg>(0,0,egg());
        }
};