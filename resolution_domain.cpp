#include "resolution_domain.hpp"
#include<iostream>

using namespace std;
using namespace meep;

void Cresolution_domain_2D::find_domain(int *dom, const vec&p){
  dom[0] = -1;
  dom[1] = -1;
  //the following code is quite ineffective, which is a problem
  //because it gets called very often. If you only have three domains
  //per dimension, which probably would be the standard use-case
  //you can implement this a bit quicker. Feel free to do so.
  //Also, I do it explicitly for two dimensions. This is ugly and
  //the whole thing can be generalized.
  
  //are you the first domain?
  if(p.x() < domain_bounds_x[0]) dom[0]=0;
  //are you the last?
  //bear in mind: dom_bounds have 1 less element than stretching!
  else if(p.x() >= domain_bounds_x[num_dom_x-2]) dom[0]=num_dom_x-1;
  //which of the inbetween are you?
  else{
    for(int i=0 ; i<(num_dom_x-2) ; ++i){
      if(p.x() >= domain_bounds_x[i] && p.x() < domain_bounds_x[i+1]){
	dom[0]=i+1;
	break;
      }
    }
  }
  //rather unnecessary failsafe:
  if(dom[0]<0)
    cout << "ERROR: did not find a domain for x= " << p.x() << endl;

  //same for the other dimension
  if(p.y() < domain_bounds_y[0]) dom[1]=0;
  else if(p.y() >= domain_bounds_y[num_dom_y-2]) dom[1]=num_dom_y-1;
  else{
    for(int i=0 ; i<(num_dom_y-2) ; ++i){
      if(p.y() >= domain_bounds_y[i] && p.y() < domain_bounds_y[i+1]){
	dom[1]=i+1;
	break;
      }
    }
  }
  if(dom[1]<0)
    cout << "ERROR: did not find a domain for y= " << p.y() << endl;
}

double Cresolution_domain_2D::chi1p1_multiplicator(const vec &p, const direction dir){
  int dom[2];
  find_domain(dom,p);
  if(dir==X)
    return stretching_x[dom[0]]/stretching_y[dom[1]];
  if(dir==Y)
    return stretching_y[dom[1]]/stretching_x[dom[0]];
  if(dir==Z)
    return 1./(stretching_x[dom[0]]*stretching_y[dom[1]]);
  
  cout << "ERROR: unexpected component direction when calculating stretch factor:"
       << dir << endl;
  return 1.0;
}
double Cresolution_domain_2D::J_multiplicator(const vec &p, const direction dir){
  int dom[2];
  //the next line is really weird, because when meep places the sources
  //it somehow uses a centered coordinate system. It doesn't do that
  //anywhere else
  find_domain(dom,p+center);
  
  if(dir==X)
    return 1./stretching_y[dom[1]];
  
  if(dir==Y)
    return 1./stretching_x[dom[0]];
  
  if(dir==Z)
    return 1./(stretching_x[dom[0]]*stretching_y[dom[1]]);
  
  cout << "ERROR: unexpected component direction when calculating stretch factor:"
       << dir << endl;
  return 1.0;
}

void Cresolution_domain_3D::find_domains(int *domains, const vec &p){
  domains[0] = domains[1] = domains[2] = -1;
  //X-direction:
  //first domain?
  if(p.x()<domain_bounds_x[0]) domains[0]=0;
  //last domain? (there is one less domain bound than there are domains)
  else if(p.x()>=domain_bounds_x[num_dom_x-2]) domains[0]=num_dom_x-1;
  //anything in between?
  else{
    for (int i=1;i<num_dom_x-1;++i){
      if(p.x() >= domain_bounds_x[i-1] && p.x() < domain_bounds_x[i]){
	domains[0]=i;
	break;
      }
    }
  }

  if(p.y()<domain_bounds_y[0]) domains[1]=0;
  //last domain? (there is one less domain bound than there are domains)
  else if(p.y()>=domain_bounds_y[num_dom_y-2]) domains[1]=num_dom_y-1;
  //anything in between?
  else{
    for (int i=1;i<num_dom_y-1;++i){
      if(p.y() >= domain_bounds_y[i-1] && p.y() < domain_bounds_y[i]){
	domains[1]=i;
	break;
      }
    }
  }

  if(p.z()<domain_bounds_z[0]) domains[2]=0;
  //last domain? (there is one less domain bound than there are domains)
  else if(p.z()>=domain_bounds_z[num_dom_z-2]) domains[2]=num_dom_z-1;
  //anything in between?
  else{
    for (int i=1;i<num_dom_z-1;++i){
      if(p.z() >= domain_bounds_z[i-1] && p.z() < domain_bounds_z[i]){
	domains[2]=i;
	break;
      }
    }
  }
  for(int i=0;i<3;++i){
    if(domains[i]<0)  cout << "ERROR! domain not found" << endl;
  }
  return;
}


double Cresolution_domain_3D::chi1p1_multiplicator(const vec &p, const direction dir){
  int cur_dom[3]; //current domain
  find_domains(cur_dom,p);
  if(dir==X)
    return stretching_x[cur_dom[0]]/(stretching_y[cur_dom[1]]*stretching_z[cur_dom[2]]);
  if(dir==Y)
    return stretching_y[cur_dom[1]]/(stretching_x[cur_dom[0]]*stretching_z[cur_dom[2]]);
  if(dir==Z)
    return stretching_z[cur_dom[2]]/(stretching_y[cur_dom[1]]*stretching_x[cur_dom[0]]);
  cout << "ERROR: unexpected component direction when calculating stretch factor:"
       << dir << endl;
  return 1.0;
}


double Cresolution_domain_3D::J_multiplicator(const vec &p, const direction dir){
  return 1.0;
  int cur_dom[3]; //current domain
  find_domains(cur_dom,p+center); //weird shift for meeps coordinate system for sources
  if(dir==X)
    return 1./(stretching_y[cur_dom[1]]*stretching_z[cur_dom[2]]);
  if(dir==Y)
    return 1./(stretching_x[cur_dom[0]]*stretching_z[cur_dom[2]]);
  if(dir==Z)
    return 1./(stretching_y[cur_dom[1]]*stretching_x[cur_dom[0]]);
  cout << "ERROR: unexpected component direction when calculating stretch factor:"
       << dir << endl;
  return 1.0;
}
