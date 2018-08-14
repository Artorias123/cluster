#include "point.h"
#include <iostream>
using namespace std;
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
//point
point::point(){
    dim=0;
}
point::point(unsigned d){
	dim=d;
	x=(double*)malloc(sizeof(double)*d);
}
point::point(const point &y){
	dim=y.dim;
	x=(double*)malloc(sizeof(double)*dim);
	for(unsigned i=0;i<dim;i++) x[i]=y[i];
}
point::point(point &&y){
	dim=y.dim;
	x=y.x;
	y.x=NULL;
}
point::point(double *y,unsigned d){
	dim=d;
	x=(double*)malloc(sizeof(double)*d);
	for(unsigned i=0;i<dim;i++) x[i]=y[i];
}
point::point(const std::vector<double> &y){
	dim=y.size();
	x=(double*)malloc(sizeof(double)*dim);
	for(unsigned i=0;i<dim;i++) x[i]=y[i];
}
point::~point(){
	if(x) {free(x);x=NULL;}
}
double& point::operator[](const unsigned &i){
	if(i<dim) return x[i];
	else{
		std::cout<<"over max"<<std::endl;
		return x[0];
	}
}
const double& point::operator[](const unsigned &i) const{
	if(i<dim) return x[i];
	else{
		std::cout<<"over max"<<std::endl;
		return x[0];
	}
}
point& point::operator=(const point &y){
	if(this!=&y){
		if(dim!=y.dim){
			dim=y.dim;
			if(x) free(x);
			x=(double*)malloc(sizeof(double)*dim);
		}
		for(unsigned i=0;i<dim;i++) x[i]=y[i];
	}
	return *this;
}
void point::show(){
	for(unsigned i=0;i<dim;i++) std::cout<<x[i]<<'\t';
	std::cout<<std::endl;
}
void point::prt(std::ofstream &fout){
	for(unsigned i=0;i<dim;i++) fout<<x[i]<<'\t';
	fout<<std::endl;
}
unsigned& point::get_dim(){
    return dim;
}
const unsigned& point::get_dim() const{
    return dim;
}
void point::zero(){
    for(unsigned i=0;i<dim;i++) x[i]=0;
}
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
//centor point
cpoint::cpoint():point(){num=0;}
cpoint::cpoint(unsigned d):point(d){for(unsigned i=0;i<dim;i++);}
cpoint::cpoint(const point &y):point(y){}
cpoint::cpoint(point &&y):point(y){}
cpoint::cpoint(double *y, unsigned d):point(y,d){}
cpoint::cpoint(const std::vector<double> &y):point(y){}
cpoint::cpoint(double *y, unsigned d, int n):point(y,d){num=n;}
cpoint::cpoint(const std::vector<double> &y, int n):point(y){num=n;}
cpoint::cpoint(const point &y, int n):point(y){num=n;}
cpoint::cpoint(point &&y, int n):point(y){num=n;}
cpoint::cpoint(const cpoint &y){
    dim=y.dim;
    num=y.num;
    x=(double*)malloc(sizeof(double)*dim);
    for(unsigned i=0;i<dim;i++) x[i]=y[i];
}
cpoint::cpoint(cpoint &&y){
    dim=y.dim;
    num=y.num;
    x=y.x;
    y.x=NULL;
};
double& cpoint::operator[](const unsigned &i){
    if(i<dim) return x[i];
	else{
		std::cout<<"over max"<<std::endl;
		return x[0];
	}
}
const double& cpoint::operator[](const unsigned &i) const{
    if(i<dim) return x[i];
	else{
		std::cout<<"over max"<<std::endl;
		return x[0];
	}
}
unsigned& cpoint::get_num(){return num;}
const unsigned& cpoint::get_num() const{return num;}
cpoint& cpoint::operator=(const cpoint &y){
    if(this!=&y){
        num=y.num;
		if(dim!=y.dim){
			dim=y.dim;
			if(x) free(x);
			x=(double*)malloc(sizeof(double)*dim);
		}
		for(unsigned i=0;i<dim;i++) x[i]=y[i];
	}
	return *this;
}
void cpoint::show(){
	for(unsigned i=0;i<dim;i++) std::cout<<x[i]<<'\t';
	std::cout<<num<<std::endl;
}
void cpoint::prt(std::ofstream &fout){
	for(unsigned i=0;i<dim;i++) fout<<x[i]<<'\t';
	fout<<num<<std::endl;
}
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
//K-Means
K_point::K_point():point(){}
K_point::K_point(unsigned d):point(d){}
K_point::K_point(const point &y):point(y){}
K_point::K_point(point &&y):point(y){}
K_point::K_point(double *y, unsigned d):point(y,d){}
K_point::K_point(const std::vector<double> &y):point(y){}
K_point::K_point(double *y, unsigned d, int k):point(y,d){kind=k;}
K_point::K_point(const std::vector<double> &y, int k):point(y){kind=k;}
K_point::K_point(const point &y, int k):point(y){kind=k;}
K_point::K_point(point &&y, int k):point(y){kind=k;}
K_point::K_point(const K_point &y){
    dim=y.dim;
    kind=y.kind;
    x=(double*)malloc(sizeof(double)*dim);
    for(unsigned i=0;i<dim;i++) x[i]=y[i];
}
K_point::K_point(K_point &&y){
    dim=y.dim;
    kind=y.kind;
    x=y.x;
    y.x=NULL;
}
K_point& K_point::operator=(const K_point &y){
	if(this!=&y){
        kind=y.kind;
		if(dim!=y.dim){
			dim=y.dim;
			if(x) free(x);
			x=(double*)malloc(sizeof(double)*dim);
		}
		for(unsigned i=0;i<dim;i++) x[i]=y[i];
	}
	return *this;
}
double& K_point::operator[](const unsigned &i){
	if(i<dim) return x[i];
	else{
		std::cout<<"over max"<<std::endl;
		return x[0];
	}
}
const double& K_point::operator[](const unsigned &i) const{
	if(i<dim) return x[i];
	else{
		std::cout<<"over max"<<std::endl;
		return x[0];
	}
}
void K_point::show(){
	for(unsigned i=0;i<dim;i++) std::cout<<x[i]<<'\t';
	std::cout<<kind<<std::endl;
}
void K_point::prt(std::ofstream &fout){
	for(unsigned i=0;i<dim;i++) fout<<x[i]<<'\t';
	fout<<kind<<std::endl;
}
unsigned& K_point::get_kind(){
    return kind;
}
const unsigned& K_point::get_kind() const{
    return kind;
}
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
//Fuzzy C-Means & Gaussian Mixed Model
FG_point::FG_point():point(){kind_n=0;}
FG_point::FG_point(unsigned d, int sk):point(d){kind_n=sk;u=(double*)malloc(sizeof(double)*sk);}
FG_point::FG_point(const point &y, int sk):point(y){kind_n=sk;u=(double*)malloc(sizeof(double)*sk);}
FG_point::FG_point(point &&y, int sk):point(y){kind_n=sk;u=(double*)malloc(sizeof(double)*sk);}
FG_point::FG_point(double *y, unsigned d, int sk):point(y,d){kind_n=sk;u=(double*)malloc(sizeof(double)*sk);}
FG_point::FG_point(const std::vector<double> &y,int sk):point(y){kind_n=sk;u=(double*)malloc(sizeof(double)*sk);}
FG_point::FG_point(double *y, unsigned d, int k, int sk):point(y,d){kind_n=sk;u=(double*)malloc(sizeof(double)*sk);kind=k;}
FG_point::FG_point(const std::vector<double> &y, int k, int sk):point(y){kind_n=sk;u=(double*)malloc(sizeof(double)*sk);kind=k;}
FG_point::FG_point(const FG_point &y){
    dim=y.dim;
    kind=y.kind;
    kind_n=y.kind_n;
    x=(double*)malloc(sizeof(double)*dim);
    u=(double*)malloc(sizeof(double)*kind_n);
    for(unsigned i=0;i<dim;i++) x[i]=y[i];
    for(unsigned i=0;i<kind_n;i++) u[i]=y.u[i];
    }
FG_point::FG_point(FG_point &&y){
    dim=y.dim;
    kind=y.kind;
    kind_n=kind_n;
    x=y.x;
    y.x=NULL;
    u=y.u;
    y.u=NULL;
}
FG_point::~FG_point(){if(!u) free(u);u=NULL;}
void FG_point::cal_kind(){
    double max=u[0];
    unsigned maxn=0;
    for(unsigned i=1;i<kind_n;i++){
        if(max<u[i]) {
            maxn=i;
            max=u[i];
        }
    }
    kind=maxn;
}
double& FG_point::operator[](const unsigned &i){
	if(i<dim) return x[i];
	else{
		std::cout<<"over max"<<std::endl;
		return x[0];
	}
}
const double& FG_point::operator[](const unsigned &i) const{
	if(i<dim) return x[i];
	else{
		std::cout<<"over max"<<std::endl;
		return x[0];
	}
}
void FG_point::show(){
	for(unsigned i=0;i<dim;i++) std::cout<<x[i]<<'\t';
	std::cout<<kind<<std::endl;
}
void FG_point::prt(std::ofstream &fout){
	for(unsigned i=0;i<dim;i++) fout<<x[i]<<'\t';
	fout<<kind<<std::endl;
}
void FG_point::show_u(){
	for(unsigned i=0;i<kind_n;i++) std::cout<<u[i]<<'\t';
	std::cout<<std::endl;
}
void FG_point::prt_u(std::ofstream &fout){
	for(unsigned i=0;i<kind_n;i++) fout<<u[i]<<'\t';
	fout<<std::endl;
}
FG_point& FG_point::operator=(const FG_point &y){
	if(this!=&y){
        kind=y.kind;
		if(dim!=y.dim){
			dim=y.dim;
			if(x) free(x);
			x=(double*)malloc(sizeof(double)*dim);
		}
		if(kind_n!=y.kind_n){
			kind_n=y.kind_n;
			if(u) free(u);
			u=(double*)malloc(sizeof(double)*kind_n);
		}
		for(unsigned i=0;i<dim;i++) x[i]=y[i];
		for(unsigned i=0;i<kind_n;i++) u[i]=y.u[i];
	}
	return *this;
}
unsigned& FG_point::get_kindn(){
	return kind_n;
}
const unsigned& FG_point::get_kindn() const{
	return kind_n;
}
unsigned& FG_point::get_kind(){
	return kind;
}
const unsigned& FG_point::get_kind() const{
	return kind;
}
double& FG_point::getu(const unsigned &i){
	if(i<kind_n) return u[i];
	else {
		std::cout<<"over u limit"<<std::endl;
		return u[0];
	} 
}
const double& FG_point::getu(const unsigned &i) const{
	if(i<kind_n) return u[i];
	else {
		std::cout<<"over u limit"<<std::endl;
		return u[0];
	} 
}
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
//Fast Search and Find of Density Peaks
DP_point::DP_point():point(){}
DP_point::DP_point(unsigned d):point(d){}
DP_point::DP_point(const point &y):point(y){}
DP_point::DP_point(point &&y):point(y){}
DP_point::DP_point(double *y, unsigned d):point(y,d){}
DP_point::DP_point(const std::vector<double> &y):point(y){}
DP_point::DP_point(const point &y, int r, double d):point(y){rho=r;delta=d;}
DP_point::DP_point(point &&y, int r, double d):point(y){rho=r;delta=d;}
DP_point::DP_point(const K_point &y){
	dim=y.get_dim();
	x=(double*)malloc(sizeof(double)*dim);
    for(unsigned i=0;i<dim;i++) x[i]=y[i];
}
DP_point::DP_point(const FG_point &y){
	dim=y.get_dim();
	x=(double*)malloc(sizeof(double)*dim);
    for(unsigned i=0;i<dim;i++) x[i]=y[i];
}
DP_point::DP_point(const DP_point &y){
    dim=y.dim;
    rho=y.rho;
    delta=y.delta;
    x=(double*)malloc(sizeof(double)*dim);
    for(unsigned i=0;i<dim;i++) x[i]=y[i];
}
DP_point::DP_point(DP_point &&y){
    dim=y.dim;
    rho=y.rho;
    delta=y.delta;
    x=y.x;
    y.x=NULL;
}
DP_point& DP_point::operator=(const DP_point &y){
	if(this!=&y){
        rho=y.rho;
        delta=y.delta;
		if(dim!=y.dim){
			dim=y.dim;
			if(x) free(x);
			x=(double*)malloc(sizeof(double)*dim);
		}
		if(!x) x=(double*)malloc(sizeof(double)*dim);
		for(unsigned i=0;i<dim;i++) x[i]=y[i];
	}
	return *this;
}
DP_point& DP_point::operator=(const K_point &y){
	if(dim!=y.get_dim()){
		dim=y.get_dim();
		if(x) free(x);
		x=(double*)malloc(sizeof(double)*dim);
	}
	if(!x) x=(double*)malloc(sizeof(double)*dim);
	for(unsigned i=0;i<dim;i++) x[i]=y[i];
	return *this;
}
DP_point& DP_point::operator=(const FG_point &y){
	if(dim!=y.get_dim()){
		dim=y.get_dim();
		if(x) free(x);
		x=(double*)malloc(sizeof(double)*dim);
	}
	if(!x) x=(double*)malloc(sizeof(double)*dim);
	for(unsigned i=0;i<dim;i++) x[i]=y[i];
	return *this;
}
double& DP_point::operator[](const unsigned &i){
	if(i<dim) return x[i];
	else{
		std::cout<<"over max"<<std::endl;
		return x[0];
	}
}
const double& DP_point::operator[](const unsigned &i) const{
	if(i<dim) return x[i];
	else{
		std::cout<<"over max"<<std::endl;
		return x[0];
	}
}
void DP_point::show(){
	for(unsigned i=0;i<dim;i++) std::cout<<x[i]<<'\t';
	std::cout<<rho<<'\t'<<delta<<std::endl;
}
void DP_point::prt(std::ofstream &fout){
	for(unsigned i=0;i<dim;i++) fout<<x[i]<<'\t';
	fout<<rho<<'\t'<<delta<<std::endl;
}
unsigned& DP_point::get_rho(){
    return rho;
}
double& DP_point::get_delta(){
    return delta;
}
const unsigned& DP_point::get_rho() const{
    return rho;
}
const double& DP_point::get_delta() const{
    return delta;
}
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
CGMM_point::CGMM_point():cpoint(){}
CGMM_point::CGMM_point(unsigned d):cpoint(d){
	sigma=(double*)malloc(sizeof(double)*d*d);
	for(unsigned i=0;i<d*d;i++) sigma[i]=0;
	for(unsigned i=0;i<d;i++) sigma[i+i*d]=250000;
};
CGMM_point::CGMM_point(const CGMM_point &y){
	dim=y.dim;
    num=y.num;
    x=(double*)malloc(sizeof(double)*dim);
    sigma=(double*)malloc(sizeof(double)*dim*dim);
    for(unsigned i=0;i<dim;i++) x[i]=y[i];
    for(unsigned i=0;i<dim*dim;i++) sigma[i]=y.sigma[i];
}
CGMM_point::CGMM_point(CGMM_point &&y){
	dim=y.dim;
    num=y.num;
	x=y.x;
	sigma=y.sigma;
	y.x=NULL;
	y.sigma=NULL;
}
CGMM_point::~CGMM_point(){
	if(!sigma) free(sigma);
	sigma=NULL;
}
double& CGMM_point::operator[](const unsigned &i){
	if(i<dim) return x[i];
	else {
		cout<<"over limit"<<endl;
		return x[0];
	}
}
const double& CGMM_point::operator[](const unsigned &i) const{
	if(i<dim) return x[i];
	else {
		cout<<"over limit"<<endl;
		return x[0];
	}
}
unsigned& CGMM_point::get_num(){
	return num;
}
const unsigned& CGMM_point::get_num() const{
	return num;
}
CGMM_point& CGMM_point::operator=(const CGMM_point &y){
	if(this!=&y){
		dim=y.dim;
    	num=y.num;
		if(dim!=y.dim){
			if(x) free(x);
			if(sigma) free(sigma);
			x=(double*)malloc(sizeof(double)*dim);
			sigma=(double*)malloc(sizeof(double)*dim*dim);
		}
    	for(unsigned i=0;i<dim;i++) x[i]=y[i];
    	for(unsigned i=0;i<dim*dim;i++) sigma[i]=y.sigma[i];
	}
	return *this;
}
void CGMM_point::show(){
	for(unsigned i=0;i<dim;i++) std::cout<<x[i]<<'\t';
	std::cout<<num<<std::endl;
}
void CGMM_point::prt(std::ofstream &fout){
	for(unsigned i=0;i<dim;i++) fout<<x[i]<<'\t';
	fout<<num<<std::endl;
}
double& CGMM_point::getsig(const unsigned &i,const unsigned &j){
	return sigma[i*dim+j];
}
const double& CGMM_point::getsig(const unsigned &i,const unsigned &j) const{
	return sigma[i*dim+j];
}
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
double dis(const point &a,const point &b){
    if(a.get_dim()!=b.get_dim()) return 0;
    double r=0;
    for(unsigned i=0;i<a.get_dim();i++) r+=(a[i]-b[i])*(a[i]-b[i]);
    return r;
}