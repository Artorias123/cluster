#include "calculator.h"
#include <algorithm>
#include <cmath>
using namespace std;
FSFDP::FSFDP(vector<DP_point> &t, vector<cpoint> &cen, double l){
    p = &t; 
    c = &cen; 
    lx = l*l;
}
bool lessr(const DP_point &x,const DP_point &y){
    return (x.get_rho()>y.get_rho());
}
bool lessd(const DP_point &x,const DP_point &y){
    return (x.get_delta()>y.get_delta());
}
void FSFDP::density(){
    double d,j;
	for (unsigned i = 0; i < p->size(); i++){
		p->at(i).get_rho() = 0;
		for (j = i+1; j < p->size(); j++){
			d = dis(p->at(i), p->at(j));
			if (d<lx){
				p->at(i).get_rho()++;
				p->at(j).get_rho()++;
			}
		}
	}
    sort(p->begin(), p->end(), lessr);
}
void FSFDP::delta(){
    double d;
	for (unsigned i = 0; i < p->size(); i++){
		d = dis(p->at(0), p->at(i));
		if (i==0||p->at(0).get_delta()<d){
			p->at(0).get_delta() = d;
		}
		for (unsigned j = 0; j < i; j++){
			d = dis(p->at(i), p->at(j));
			if (j==0||p->at(i).get_delta()>d) p->at(i).get_delta() = d;
		}
        p->at(i).get_delta() = pow(p->at(i).get_delta(),0.5)*double(p->at(i).get_rho());
	}
    sort(p->begin(), p->end(), lessd);
}
void FSFDP::search(){
    unsigned j,n=p->at(0).get_dim();
    for(unsigned i=0;i<c->size();i++){
        for(j=0;j<n;j++) c->at(i)[j]=p->at(i)[j];
    }
}
void FSFDP::calculator(){
    density();
	delta();
	search();
    ofstream dpfout("./result/FSFDP.dat",ios::trunc);
    for(unsigned i=0;i<c->size();i++){
        for(unsigned j=0;j<c->at(0).get_dim();j++) dpfout<<c->at(i)[j]<<'\t';
        dpfout<<endl;
    }
    dpfout.close();
}
void FSFDP::prt_result(std::ofstream &cfout){
    unsigned i=0,j=0,dim=p->at(0).get_dim();
    for(;i<c->size();i++){
        for(;j<dim;j++) cfout<<c->at(i)[j]<<'\t';
        cfout<<endl;
    }
}