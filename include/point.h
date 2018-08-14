#ifndef point_H
#define point_H
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
//point
class point{
    protected:
    double *x;
    unsigned dim;
    public:
    point();
    point(unsigned d);
    point(const point &y);
    point(point &&y);
    point(double *y,unsigned d);
    point(const std::vector<double> &y);
    ~point();
    virtual double& operator[](const unsigned &i);
    virtual const double& operator[](const unsigned &i) const;
    point& operator=(const point &y);
    virtual void show();
    virtual void prt(std::ofstream &fout);
    unsigned& get_dim();
    const unsigned& get_dim() const;
    void zero();
};
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
//centor point
class cpoint : public point{
    protected:
    unsigned num;
    public:
    cpoint();
    cpoint(unsigned d);
    cpoint(const point &y);
    cpoint(point &&y);
    cpoint(double *y, unsigned d);
    cpoint(const std::vector<double> &y);
    cpoint(double *y, unsigned d, int n);
    cpoint(const std::vector<double> &y, int n);
    cpoint(const point &y, int n);
    cpoint(point &&y, int n);
    cpoint(const cpoint &y);
    cpoint(cpoint &&y);
    virtual double& operator[](const unsigned &i);
    virtual const double& operator[](const unsigned &i) const;
    unsigned& get_num();
    const unsigned& get_num() const;
    cpoint& operator=(const cpoint &y);
    virtual void show();
    virtual void prt(std::ofstream &fout);
};
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
//K-Means
class K_point : public point{
    protected:
    unsigned kind;
    protected:
    double* get_xp();
    const double* get_xp() const;
    public:
    K_point();
    K_point(unsigned d);
    K_point(const point &y);
    K_point(point &&y);
    K_point(double *y, unsigned d);
    K_point(const std::vector<double> &y);
    K_point(double *y, unsigned d, int k);
    K_point(const std::vector<double> &y, int k);
    K_point(const point &y, int k);
    K_point(point &&y, int k);
    K_point(const K_point &y);
    K_point(K_point &&y);
    K_point& operator=(const K_point &y);
    virtual double& operator[](const unsigned &i);
    virtual const double& operator[](const unsigned &i) const;
    virtual void show();
    virtual void prt(std::ofstream &fout);
    unsigned& get_kind();
    const unsigned& get_kind() const;
};
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
//Fuzzy C-Means & Gaussian Mixed Model
class FG_point : public point{
    private:
    double *u;
    unsigned kind_n,kind;
    public:
    FG_point();
    FG_point(unsigned d, int sk);
    FG_point(const point &y, int sk);
    FG_point(point &&y, int sk);
    FG_point(double *y, unsigned d, int sk);
    FG_point(const std::vector<double> &y,int sk);
    FG_point(double *y, unsigned d, int k, int sk);
    FG_point(const std::vector<double> &y, int k, int sk);
    FG_point(const FG_point &y);
    FG_point(FG_point &&y);
    ~FG_point();
    void cal_kind();
    virtual double& operator[](const unsigned &i);
    virtual const double& operator[](const unsigned &i) const;
    virtual void show();
    virtual void prt(std::ofstream &fout);
    void show_u();
    void prt_u(std::ofstream &fout);
    FG_point& operator=(const FG_point &y);
    unsigned& get_kindn();
    const unsigned& get_kindn() const;
    unsigned& get_kind();
    const unsigned& get_kind() const;
    double& getu(const unsigned &i);
    const double& getu(const unsigned &i) const;
};
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
//Fast Search and Find of Density Peaks
class DP_point : public point{
    private:
    unsigned rho;
    double delta;
    public:
    DP_point();
    DP_point(unsigned d);
    DP_point(const point &y);
    DP_point(point &&y);
    DP_point(double *y, unsigned d);
    DP_point(const std::vector<double> &y);
    DP_point(const point &y, int r, double d);
    DP_point(point &&y, int r, double d);
    DP_point(const K_point &y);
    DP_point(const FG_point &y);
    DP_point(const DP_point &y);
    DP_point(DP_point &&y);
    DP_point& operator=(const DP_point &y);
    DP_point& operator=(const K_point &y);
    DP_point& operator=(const FG_point &y);
    virtual double& operator[](const unsigned &i);
    virtual const double& operator[](const unsigned &i) const;
    virtual void show();
    virtual void prt(std::ofstream &fout);
    unsigned& get_rho();
    double& get_delta();
    const unsigned& get_rho() const;
    const double& get_delta() const;
};
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
//Gaussian Mixed Model centor point
class CGMM_point : public cpoint{
    private:
    double *sigma;
    public:
    CGMM_point();
    CGMM_point(unsigned d);
    CGMM_point(const CGMM_point &y);
    CGMM_point(CGMM_point &&y);
    ~CGMM_point();
    virtual double& operator[](const unsigned &i);
    virtual const double& operator[](const unsigned &i) const;
    unsigned& get_num();
    const unsigned& get_num() const;
    CGMM_point& operator=(const CGMM_point &y);
    virtual void show();
    virtual void prt(std::ofstream &fout);
    double& getsig(const unsigned &i,const unsigned &j);
    const double& getsig(const unsigned &i,const unsigned &j) const;
};
//-------1---------2---------3---------4---------5---------6---------7---------8---------9---------
double dis(const point &a,const point &b);
void c2CG(std::vector<CGMM_point> &c,const std::vector<cpoint> &tc,const unsigned &np);
#endif
