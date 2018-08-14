#include "calculator.h"
#include <cmath>
using namespace std;
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
FCMeans::FCMeans(vector<FG_point> &t, vector<cpoint> &cen,double im, unsigned ni):method(ni){
    p = &t; 
    c = &cen; 
    oldc=cen;
    m=im;
}
void FCMeans::init(){
    unsigned tr,j,dim=c->at(0).get_dim(),n=p->size();
    vector<unsigned> tmp(n);
    for(j=0;j<n;j++) tmp[j]=j;
    for(unsigned i=0;i<c->size();i++){
        tr=rand()%tmp.size();
        for(j=0;j<dim;j++) c->at(i)[j]=p->at(tmp[tr])[j];
        tmp.erase(tmp.begin()+tr);
    }
}
void FCMeans::update_p(){
    oldc=*c;
    unsigned i,j,k,n=c->size();
    double tmp,td;
    for (i = 0; i < n; i++){
		for (j = 0; j < p->size(); j++){
			tmp = 0;
			for (k = 0; k < n; k++) {
                td=dis(p->at(j), c->at(k));
                if(td<pow(0.1,12)) td=pow(0.1,12);
                tmp += pow(dis(p->at(j), c->at(i)) / td, 2 / (m - 1));
            }
            if(tmp<pow(0.1,12)) p->at(j).getu(i)=1;
			else p->at(j).getu(i) = 1 / tmp;
		}
	}
}
void FCMeans::update_cen(){
    unsigned i,j,k,n=c->at(0).get_dim();
    double tmp;
    for (i = 0; i < c->size(); i++) c->at(i).zero();
    for (i = 0; i < c->size(); i++) {
		tmp=0;
        for (j = 0; j < p->size(); j++){
			for(k=0;k<n;k++) c->at(i)[k] += pow(p->at(j).getu(i), m) * p->at(j)[k];
            tmp += pow(p->at(j).getu(i), m);
		}
        for(k=0;k<n;k++) c->at(i)[k] /=tmp;
	}
}
double FCMeans::cal_J(){
    if(step==0) return 100;
    double J=0;
    for(unsigned i=0;i<c->size();i++) J+=dis(c->at(i),oldc[i]);
    return J;
}
void FCMeans::prt_result(std::ofstream &pfout,std::ofstream &cfout){
    classify();
    unsigned i,j,dim=p->at(0).get_dim();
    for(i=0;i<p->size();i++){
        for(j=0;j<dim;j++) pfout<<p->at(i)[j]<<'\t';
        pfout<<p->at(i).get_kind()<<endl;
    }
    for(i=0;i<c->size();i++){
        for(j=0;j<dim;j++) cfout<<c->at(i)[j]<<'\t';
        cfout<<c->at(i).get_num()<<endl;
    }
}
void FCMeans::classify(){
    unsigned i=0,j,n=c->size();
    double max;
    for(;i<n;i++) c->at(i).get_num()=0;
    for(i=0;i<p->size();i++){
        for(j=0;j<n;j++){
            if(j==0||max<p->at(i).getu(j)) {
                max=p->at(i).getu(j);
                p->at(i).get_kind()=j;
            }
        }
        c->at(p->at(i).get_kind()).get_num()++;
    }
}