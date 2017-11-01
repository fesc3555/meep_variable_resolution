#include "material.hpp"

vec Ceps::normal_vector_anis(field_type ft, const volume &v,const direction dir){
  vec gradient(zero_vec(v.dim));
  vec p(v.center());
  double R = v.diameter();
  for (int i = 0; i < num_sphere_quad[number_of_directions(v.dim)-1]; ++i) {
    double weight;
    vec pt = sphere_pt(p, R, i, weight);
    //probably only a good formula for diagonal chi1p1:
    gradient += (pt - p) * (weight * res_dom.chi1p1_multiplicator(pt,dir)*
			    chi1p1(ft,pt));
  }
  return gradient; 
}
  
void Ceps::eff_chi1inv_row(component c, double chi1inv_row[3], const volume &v,
			   double tol, int maxeval){
  field_type ft = type(c);
  vec gradient(normal_vector_anis(ft, v,component_direction(c)));
  if (!maxeval || abs(gradient) < 1e-8 ){
  trivial:
    chi1inv_row[0] = chi1inv_row[1] = chi1inv_row[2] = 0.;
    //chi1pi = J*chi1p1*J'/det(J) boils down to just a factor for diagonal eps
    chi1inv_row[component_direction(c)]=
      1/(res_dom.chi1p1_multiplicator(v.center(),component_direction(c)) *
	 chi1p1(ft,v.center()));
    return;
  }
  double meps=1, minveps=1;
  vec d = v.get_max_corner() - v.get_min_corner();
  int ms = 5;
  double old_meps=0, old_minveps=0;
  int iter = 0;
  while ((fabs(meps-old_meps) > tol*old_meps) && (fabs(minveps-old_minveps) > tol*old_minveps)) {
    old_meps=meps; old_minveps=minveps;
    meps = minveps = 0;
    for (int j=0; j < ms; j++)
      for (int i=0; i < ms; i++) {
	vec p=v.get_min_corner() + vec(i*d.x()/ms, j*d.y()/ms);
	double ep = res_dom.chi1p1_multiplicator(p,component_direction(c))*chi1p1(ft,p);
	if (ep < 0) {goto trivial;}
	meps += ep; minveps += 1/ep;
      }
    meps /= ms*ms;
    minveps /= ms*ms;
    ms *= 2;
    if (maxeval && (iter += ms*ms) >= maxeval) break;
  }
  double n[3] = {0,0,0};
  double nabsinv = 1.0/abs(gradient);
  LOOP_OVER_DIRECTIONS(gradient.dim, k)
    n[k%3] = gradient.in_direction(k) * nabsinv;
  
  /* get rownum'th row of effective tensor
     P * minveps + (I-P) * 1/meps = P * (minveps-1/meps) + I * 1/meps
     where I is the identity and P is the projection matrix
     P_{ij} = n[i] * n[j]. */
  int rownum = component_direction(c) % 3;
  for (int i=0; i<3; ++i)
    chi1inv_row[i] = n[rownum] * n[i] * (minveps - 1/meps);
  chi1inv_row[rownum] += 1/meps;
}

