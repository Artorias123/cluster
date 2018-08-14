#include "calculator.h"
using namespace std;
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
KMeans::KMeans(vector<K_point> &t, vector<cpoint> &cen, unsigned ni):method(ni){
    p = &t; 
    c = &cen; 
    oldc=cen;
}
void KMeans::init(){
    unsigned tr,j,dim=c->at(0).get_dim(),n=p->size();
    vector<unsigned> tmp(n);
    for(j=0;j<n;j++) tmp[j]=j;
    for(unsigned i=0;i<c->size();i++){
        tr=rand()%tmp.size();
        for(j=0;j<dim;j++) c->at(i)[j]=p->at(tmp[tr])[j];
        tmp.erase(tmp.begin()+tr);
    }
}
void KMeans::update_cen(){
    oldc=*c;
    unsigned i,j,n=p->at(0).get_dim();
    for(i=0;i<c->size();i++) c->at(i).zero();
    for(i=0;i<p->size();i++){
        for(j=0;j<n;j++) c->at(p->at(i).get_kind())[j]+=p->at(i)[j];
    }
    for(i=0;i<c->size();i++) {
        if(c->at(i).get_num()!=0) for(j=0;j<n;j++) c->at(i)[j]/=c->at(i).get_num();
        else for(j=0;j<n;j++) c->at(i)[j]=p->at(rand()%p->size())[j];
    }
}
void KMeans::update_p(){
    unsigned i,j,n=c->size(),minn=0;
    double min,tmp;
    for(j=0;j<n;j++) c->at(j).get_num()=0;
    for(i=0;i<p->size();i++){
        for(j=0;j<n;j++){
            tmp=dis(p->at(i),c->at(j));
            if(j==0||min>tmp){
                min=tmp;
                minn=j;
            }
        }
        c->at(minn).get_num()++;
        p->at(i).get_kind()=minn;
    }
}
double KMeans::cal_J(){
    if(step==0) return 100;
    double J=0;
    for(unsigned i=0;i<c->size();i++) J+=dis(c->at(i),oldc[i]);
    return J;
}
void KMeans::prt_result(std::ofstream &pfout,std::ofstream &cfout){
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