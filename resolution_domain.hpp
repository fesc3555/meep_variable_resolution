#ifndef resdom_fesc
#define resdom_fesc

#include<meep.hpp>
#include<vector>

using namespace meep;
using namespace std;

class Cresolution_domain_2D
{
public:
  Cresolution_domain_2D(vector<double> _domain_bounds_x, vector<double> _domain_bounds_y,
			vector<double> _stretching_x, vector<double> _stretching_y,
			const vec _center):
    stretching_x(_stretching_x), stretching_y(_stretching_y),
    domain_bounds_x(_domain_bounds_x),domain_bounds_y(_domain_bounds_y),
    center(_center){
    num_dom_x=stretching_x.size();
    num_dom_y=stretching_y.size();
  }
  virtual void find_domain(int *dom,const vec &p);
  double chi1p1_multiplicator(const vec &p, const direction dir);
  double J_multiplicator(const vec &p_dmy, const direction dir);
private:
  vector<double> stretching_x, stretching_y, domain_bounds_x, domain_bounds_y;
  int num_dom_x, num_dom_y;
  vec center;
};

class Cresolution_domain_3D{
public:
  Cresolution_domain_3D(vector<double> _domain_bounds_x,
		     vector<double> _domain_bounds_y,
		     vector<double> _domain_bounds_z,
		     vector<double> _stretching_x,
		     vector<double> _stretching_y,
		     vector<double> _stretching_z,
		     vec _center):
    stretching_x(_stretching_x),
    stretching_y(_stretching_y),
    stretching_z(_stretching_z),
    domain_bounds_x(_domain_bounds_x),
    domain_bounds_y(_domain_bounds_y),
    domain_bounds_z(_domain_bounds_z),
    center(_center)
  {
    num_dom_x=stretching_x.size();
    num_dom_y=stretching_y.size();
    num_dom_z=stretching_z.size();
  }
  void find_domains(int *dom, const vec &p);
  double chi1p1_multiplicator(const vec &p, const direction dir);
  double J_multiplicator(const vec &p, const direction dir);
private:
  vector<double> stretching_x, stretching_y, stretching_z,
    domain_bounds_x, domain_bounds_y, domain_bounds_z;
  int num_dom_x, num_dom_y, num_dom_z;
  vec center;
};

#endif
