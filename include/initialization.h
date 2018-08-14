#ifndef initialization_H
#define initialization_H
#include <string>
#include <fstream>
#include <vector>
#include "point.h"
enum meth {K_Means=0,FC_Means=1,GM_M=2};
class readset{
    public:
    meth m;
    std::string path,name;
    double presion,judge_dis;
    unsigned cluster_num,niter,dim,point_num;
    bool fsfdp;
    std::ifstream fin;
    std::ofstream pfout,cfout;
    readset();
    ~readset();
};
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
void readdata(std::ifstream &fin,std::vector<K_point> &p,unsigned end);
void readdata(std::ifstream &fin,std::vector<FG_point> &p,unsigned end);
void init_centor(std::vector<cpoint> &c,const std::vector<double> &lim);
void init_centor(std::vector<CGMM_point> &c,const std::vector<double> &lim,const unsigned &np);
#endif