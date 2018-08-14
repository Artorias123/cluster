#include "initialization.h"
#include "point.h"
using namespace std;
extern enum meth m;
readset::readset(){
    fin.open("set.txt");
    pfout.open("./result/point.dat",ios::trunc);
    cfout.open("./result/centor.dat",ios::trunc);
    string line,oldline;
    dim=0;
    while(getline(fin,line)){
        if(oldline=="FILENAME:") name=line;
        if(oldline=="PATH:")     path=line;
        if(oldline=="METHOD:"){
            switch(stoi(line)){
                case 0:{m=K_Means;break;}
                case 1:{m=FC_Means;break;}
                case 2:{m=GM_M;break;}
                default:{cout<<"method num error,have used K-Means";m=K_Means;}
            }
        }
        if(oldline=="DIMENSION:") dim=stoi(line);
        if(oldline=="PRESION:") presion=stod(line);
        if(oldline=="FSFDP:") fsfdp=stoi(line);
        if(oldline=="CLUSTER NUMBER:") cluster_num=stoi(line);
        if(oldline=="ITERATION NUMBER:") niter=stoi(line);
        if(oldline=="POINT NUMBER:") point_num=stoi(line);
        if(oldline=="FSFDP JUDGE DISTANCE:") judge_dis=stod(line);
        oldline=line;
    }
    fin.close();
    path=path+"/"+name;
    fin.open(path.data());
}
readset::~readset(){
    if(!fin) fin.close();
}
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
void readdata(ifstream &fin,vector<K_point> &p,unsigned end){
    unsigned n=p[0].get_dim();
    for(unsigned i=0;i<end-1;i++){
        for(unsigned j=0;j<n;j++) {
            fin>>p[i][j];
        }
        p.push_back(K_point(n));
    }
    for(unsigned j=0;j<n;j++) {
            fin>>p[p.size()-1][j];
        }
}
void readdata(ifstream &fin,vector<FG_point> &p,unsigned end){
    unsigned n=p[0].get_dim();
    for(unsigned i=0;i<end-1;i++){
        for(unsigned j=0;j<n;j++) {
            fin>>p[i][j];
        }
        FG_point tmp(n,p[0].get_kindn());
        p.push_back(tmp);
        p[i+1]=tmp;
    }
    for(unsigned j=0;j<n;j++) {
            fin>>p[p.size()-1][j];
        }
}
void init_centor(vector<cpoint> &c,const vector<double> &lim){
    unsigned i,j,n=c[0].get_dim();
    for(i=0;i<c.size();i++){
        for(j=0;j<n;j++) c[i][j]=rand()*lim[j]/RAND_MAX;
    }
}
void init_centor(vector<CGMM_point> &c,const vector<double> &lim,const unsigned &np){
    for(unsigned i=0;i<c.size();i++){
        
    }
}