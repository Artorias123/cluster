#include "calculator.h"
#include "initialization.h"
#include "point.h"
#include <cstdlib>
#include <ctime>
using namespace std;
int main(){
    srand(clock());
    readset setup;
    switch(setup.m){
        case 0:{
            K_point x(setup.dim);
            cpoint y(setup.dim);
            vector<K_point> p(1,x);
            vector<cpoint> c(setup.cluster_num,y);
            p.reserve(setup.point_num);
            readdata(setup.fin,p,setup.point_num);
            if(setup.fsfdp) {
                DP_point z(x);
                vector<DP_point> tp(5000,z);
                for(unsigned i=0;i<setup.point_num-1;i++) {
                    tp[i]=p[i];
                }
                FSFDP tmp(tp,c,setup.judge_dis);
                tmp.calculator();
            }
            KMeans cal(p,c,setup.niter);
            cal.iterate(!setup.fsfdp);
            cal.prt_result(setup.pfout,setup.cfout);
            break;
        }
        case 1:{
            FG_point x(setup.dim,setup.cluster_num);
            cpoint y(setup.dim);
            vector<FG_point> p(1,x);
            vector<cpoint> c(setup.cluster_num,y);
            p.reserve(setup.point_num);
            readdata(setup.fin,p,setup.point_num);
            if(setup.fsfdp) {
                DP_point z(x);
                vector<DP_point> tp(5000,z);
                for(unsigned i=0;i<setup.point_num-1;i++) {
                    tp[i]=p[i];
                }
                FSFDP tmp(tp,c,setup.judge_dis);
                tmp.calculator();
            }
            FCMeans cal(p,c,2,setup.niter);
            cal.iterate(!setup.fsfdp);
            cal.prt_result(setup.pfout,setup.cfout);
            break;
        }
        case 2:{
            FG_point x(setup.dim,setup.cluster_num);
            CGMM_point y(setup.dim);
            vector<FG_point> p(1,x);
            vector<CGMM_point> c(setup.cluster_num,y);
            p.reserve(setup.point_num);
            readdata(setup.fin,p,setup.point_num);
            if(setup.fsfdp) {
                DP_point z(x);
                vector<cpoint> tc(setup.cluster_num,y);
                vector<DP_point> tp(5000,z);
                for(unsigned i=0;i<setup.point_num-1;i++) {
                    tp[i]=p[i];
                }
                FSFDP tmp(tp,tc,setup.judge_dis);
                tmp.calculator();
                c2CG(c,tc,setup.point_num);
            }
            GMM cal(p,c,setup.niter);
            cal.init();
            cal.iterate(!setup.fsfdp);
            cal.prt_result(setup.pfout,setup.cfout);
            break;
        }
    }
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
}