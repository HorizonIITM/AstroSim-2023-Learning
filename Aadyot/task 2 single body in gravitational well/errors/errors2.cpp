
class c{
    public:
        static int d;
        // c(int d1){
        //     d=d1;
        // }

        class e1{
            int operator()() {
                return c::d;
            }
        };
};

int main(){
    c my_c;
    return 0;
}