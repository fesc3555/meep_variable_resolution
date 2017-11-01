#ifndef material_fesc
#define material_fesc
#include<meep.hpp>
#include "sphere-quad.h"
#include "resolution_domain.hpp"

using namespace std;
using namespace meep;

//this function is just copied&pasted from anisotropic_averaging.cpp in the meep src
static vec sphere_pt(const vec &cent, double R, int n, double &weight) {
     switch (cent.dim) {
         case D1:
         {
              weight = sphere_quad[0][n][3];
              vec pt(sphere_quad[0][n][2]);
              return cent + pt * R;
         }
         case D2:
         {
              weight = sphere_quad[1][n][3];
              vec pt(sphere_quad[1][n][0], sphere_quad[1][n][1]);
              return cent + pt * R;
         }
         case D3:
         {
              weight = sphere_quad[2][n][3];
              vec pt(sphere_quad[2][n][0], sphere_quad[2][n][1],
                     sphere_quad[2][n][2]);
              return cent + pt * R;
         }
         case Dcyl:
         {
              weight = sphere_quad[1][n][3];
              return cent
                + veccyl(sphere_quad[1][n][0], sphere_quad[1][n][1]) * R;
         }
         default:
           abort("unknown dimensions in sphere_pt\n");
     }
}

static double abs2D(const vec &a){return sqrt(a.x()*a.x()+a.y()*a.y());}

class Ceps : public material_function
{
public:
  Ceps(const Cresolution_domain_2D _res_dom, const double _sphere_rad,
       const double _eps_sphere, const vec &_center):
    res_dom(_res_dom),sphere_rad(_sphere_rad),
    eps_sphere(_eps_sphere),center(_center){}
  virtual bool has_mu(){return true;}
  virtual double chi1p1(field_type ft, const vec &r){
    if(ft == E_stuff && abs2D(r-center) < sphere_rad) return eps_sphere;
    return 1.0;
  }

  vec normal_vector_anis(field_type ft, const volume &v,const direction dir);
  
  virtual void eff_chi1inv_row(component c, double chi1inv_row[3], const volume &v,
			       double tol=DEFAULT_SUBPIXEL_TOL,
			       int maxeval=DEFAULT_SUBPIXEL_MAXEVAL);
private:
  Cresolution_domain_2D res_dom;
  double sphere_rad, eps_sphere;
  vec center;
};

#endif
