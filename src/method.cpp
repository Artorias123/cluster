#include "calculator.h"
using namespace std;
method::method(){
    step=0;
}
method::method(unsigned ni){
    n_iter=ni;
    step=0;
}
void method::advance(){
    update_p();
    update_cen();
}
void method::iterate(bool binit){
    if(binit) init();
    for(;step<n_iter&&cal_J()>0.00001;step++) advance();
}