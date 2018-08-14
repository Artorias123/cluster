#ifndef calculator_H
#define calculator_H
#include "point.h"
#include <fstream>
class method{
    public:
    unsigned n_iter,step;
    method();
    method(unsigned ni);
    virtual void update_cen()=0;
    virtual void update_p()=0;
    virtual double cal_J()=0;
    virtual void prt_result(std::ofstream &pfout,std::ofstream &cfout)=0;
    virtual void init()=0;
    void advance();
    void iterate(bool binit);
};
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
//K-Means
class KMeans : public method{
    private:
    std::vector<K_point> *p;
    std::vector<cpoint> *c;
    std::vector<cpoint> oldc;
    public:
    KMeans(std::vector<K_point> &t, std::vector<cpoint> &cen,unsigned ni);
    virtual void init();
    virtual void update_cen();
    virtual void update_p();
    virtual double cal_J();
    virtual void prt_result(std::ofstream &pfout,std::ofstream &cfout);
};
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
//Fuzzy C-Means
class FCMeans : public method{
    private:
    std::vector<FG_point> *p;
    std::vector<cpoint> *c;
    std::vector<cpoint> oldc;
    double m;
    public:
    FCMeans(std::vector<FG_point> &t, std::vector<cpoint> &cen,double im,unsigned ni);
    virtual void init();
    virtual void update_cen();
    virtual void update_p();
    virtual double cal_J();
    virtual void prt_result(std::ofstream &pfout,std::ofstream &cfout);
    void classify();
};
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
//Gaussian Mixed Model
class GMM : public method{
    private:
    std::vector<FG_point> *p;
    std::vector<CGMM_point> *c;
    std::vector<CGMM_point> oldc;
    double m;
    public:
    GMM(std::vector<FG_point> &t, std::vector<CGMM_point> &cen,unsigned ni);
    virtual void init();
    virtual void update_cen();
    virtual void update_p();
    virtual double cal_J();
    virtual void prt_result(std::ofstream &pfout,std::ofstream &cfout);
    void classify();
    void get_sig(unsigned k,unsigned l,double *x,double *y);
};
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
//Fast Search and Find of Density Peaks
class FSFDP{
    private:
    std::vector<DP_point> *p;
    std::vector<cpoint> *c;
    double lx;
    public:
    FSFDP(std::vector<DP_point> &t, std::vector<cpoint> &cen, double l);
    void density();
	void delta();
	void search();
	void calculator();
    void prt_result(std::ofstream &cfout);
};
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
#endif