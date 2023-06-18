class functor{
    public:
        void operator()(){}
};

template<class f>
class solver{
    public:
    solver(){}
    void solv1(f f1){}
};


class gsol{
    public:
        solver<functor> mysolver;
        gsol(){
            solver<functor> mysolver = solver<functor>();
        }
};