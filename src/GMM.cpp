#include "calculator.h"
#include <cmath>
using namespace std;
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
//matrix
void nmult(double *A,double *y,const unsigned &m,const double &k,const unsigned &n){
    for(unsigned i=(m-1)*n;i<m*n;i++) A[i]*=k;
    y[m-1]*=k;
}
void add(double *A,double *y,const unsigned &a,const unsigned &b,const double &k,const unsigned &n){
    for(unsigned i=0;i<n;i++) A[(b-1)*n+i]+=k*A[(a-1)*n+i];
    y[b-1]+=k*y[a-1];
}
void swap(double *A,double *y,const unsigned &a,const unsigned &b,const unsigned &n){
    double tmp;
    for(unsigned i=0;i<n;i++){
        tmp=A[(b-1)*n+i];
        A[(b-1)*n+i]=A[(a-1)*n+i];
        A[(a-1)*n+i]=tmp;
    }
    tmp=y[b];
    y[b]=y[a];
    y[a]=tmp;
}
double solqs(double *A,double *y,double *x,unsigned n){
    double m,d=1;
    for(unsigned i=0;i<n;i++){
        for(unsigned j=i;j<n;j++){
            if(A[i+j*n]) {
                swap(A,y,j+1,i+1,n);
                break;
            }
        }
        m=A[i+i*n];
        if(m==0) continue;
        d*=m;
        nmult(A,y,i+1,1/m,n);
        for(unsigned j=i+1;j<n;j++) add(A,y,i+1,j+1,-A[j*n+i],n);
    }
    d*=A[(n-1)*n+n-1];
    for(unsigned i=n-1;i>0;i--){
        for(int j=i-1;j>=0;j--) {
            add(A,y,i+1,j+1,-A[j*n+i],3);
        }
        x[i]=y[i];
    }
    x[0]=y[0];
    return (d);
}
double docmul(double *x,double *y,unsigned n){
    double r=0;
    for(unsigned i=0;i<n;i++) r+=x[i]*y[i];
    return r;
}
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
//Gaussian distribution
double N(const FG_point &p,const CGMM_point &c){
    double f,*x,*y,*A;
    x=(double*)malloc(sizeof(double)*p.get_dim());
    y=(double*)malloc(sizeof(double)*p.get_dim());
    A=(double*)malloc(sizeof(double)*p.get_dim()*p.get_dim());
    for(unsigned i=0;i<p.get_dim();i++){
        x[i]=p[i]-c[i];
        for(unsigned j=0;j<p.get_dim();j++) A[i*p.get_dim()+j]=c.getsig(i,j);
    }
    f=1/pow(solqs(A,x,y,p.get_dim()),0.5);
    for(unsigned i=0;i<p.get_dim();i++) x[i]=p[i]-c[i];
    f*=double(c.get_num())*exp(-docmul(x,y,p.get_dim())/2);
    free(x);
    free(y);
    free(A);
    return f;
}
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
GMM::GMM(std::vector<FG_point> &t, std::vector<CGMM_point> &cen,unsigned ni):method(ni){
    p = &t; 
    c = &cen; 
    oldc=cen;
}
void GMM::init(){
    unsigned minn;
    double min,tmp,*tmp0;
    for(unsigned i=0;i<c->size();i++) c->at(i).get_num()=0;
    for(unsigned i=0;i<p->size();i++){
        for(unsigned j=0;j<c->size();j++){
            tmp=dis(p->at(i),c->at(j));
            if(j==0||tmp<min){
                min=tmp;
                minn=j;
            }
        }
        p->at(i).get_kind()=minn;
        c->at(minn).get_num()++;
    }
    tmp0=(double*)malloc(sizeof(double)*p->at(0).get_dim());
    for(unsigned i=0;i<p->size();i++){
        for(unsigned k=0;k<p->at(0).get_kindn();k++) p->at(i).getu(k)=1;
        for(unsigned k=0;k<p->at(0).get_dim();k++) tmp0[k]=p->at(i)[k]-c->at(p->at(i).get_kind())[k];
        get_sig(p->at(i).get_kind(),i,tmp0,tmp0);
    }
    for(unsigned i=0;i<c->size();i++){
        for(unsigned j=0;j<p->at(0).get_dim();j++) {
            for(unsigned k=0;k<p->at(0).get_dim();k++) c->at(i).getsig(j,k)/=c->at(i).get_num()-1;
        }
    }
}
void GMM::update_cen(){
    unsigned j,k,dim=p->at(0).get_dim(),cn=c->size(),pn=p->size();
    double *tmp0,sum;
    for (unsigned i = 0; i < cn; i++){
		sum = 0;
		for (j = 0; j < pn; j++){
			 sum+= p->at(j).getu(i);
		}
        if(sum<2) sum=2;
		c->at(i).get_num() = round(sum);
	}
    tmp0=(double*)malloc(sizeof(double)*dim);
    for(unsigned i=0;i<cn;i++){
        c->at(i).zero();
        for(j=0;j<dim;j++) {
            for(k=0;k<dim;k++) c->at(i).getsig(j,k)=0;
        }
        for(j=0;j<pn;j++){
            for(k=0;k<dim;k++) c->at(i)[k]+=p->at(j)[k]*p->at(j).getu(i);
        }
        for(k=0;k<dim;k++) c->at(i)[k]/=c->at(i).get_num();
        for(j=0;j<pn;j++){
            for(k=0;k<dim;k++){
                tmp0[k]=p->at(j)[k]-c->at(i)[k];
            }
            get_sig(i,j,tmp0,tmp0);
        }
        for(j=0;j<dim;j++) {
            for(k=0;k<dim;k++) c->at(i).getsig(j,k)/=c->at(i).get_num()-1;
        }
    }
    free(tmp0);
}
void GMM::update_p(){
    double j,tmp,cn=c->size();
    for(unsigned i=0;i<p->size();i++){
        tmp=0;
        for(j=0;j<cn;j++) p->at(i).getu(j)=c->at(j).get_num()*N(p->at(i),c->at(j));
        for(j=0;j<cn;j++) tmp+=p->at(i).getu(j);
        for(j=0;j<cn;j++) p->at(i).getu(j)/=tmp;
    }
}
double GMM::cal_J(){
    if(step==0) return 100;
    double J=0;
    for(unsigned i=0;i<c->size();i++) J+=dis(c->at(i),oldc[i]);
    return J;
}
void GMM::prt_result(std::ofstream &pfout,std::ofstream &cfout){
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
void GMM::classify(){
    unsigned i=0,j,n=c->size();
    double max;
    for(i=0;i<p->size();i++){
        for(j=0;j<n;j++){
            if(j==0||max<p->at(i).getu(j)) {
                max=p->at(i).getu(j);
                p->at(i).get_kind()=j;
            }
        }
    }
}
void GMM::get_sig(unsigned k,unsigned l,double *x,double *y){
    unsigned dim=p->at(0).get_dim();
    for(unsigned i=0;i<dim;i++){
        for(unsigned j=i;j<dim;j++){
            c->at(k).getsig(i,j)+=x[i]*y[j]*p->at(l).getu(k);
            c->at(k).getsig(j,i)=c->at(k).getsig(i,j);
        }
    }
}
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
void c2CG(vector<CGMM_point> &c,const vector<cpoint> &tc,const unsigned &np){
    unsigned j,dim=c[0].get_dim();
    for(unsigned i=0;i<c.size();i++){
        for(j=0;j<dim;j++) c[i][j]=tc[i][j];
        c[i].get_num()=round(np/c.size());
    }
    c[c.size()-1].get_num()=np-(c.size()-1)*round(np/c.size());
}