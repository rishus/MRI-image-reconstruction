#include "adaptive_simpsons.h"
#include "adaptive_simpsons_profiled.h"

#ifdef __USE_QUANTLIB__
GaussLegendreIntegration gintegrator(64);
#endif
//AdaptiveSimponsIntegration integrator;
//AdaptiveSimponsIntegrationProfiled  pintegrator;



template<int i, class integrator>
struct I_l_n;

template<int i>
struct F_l_n;



extern Real R;
extern Real tol;




void cosntheta(vec_rel_ptr_t cptr, Real costheta, int n ) //result in cptr
{
    int k;
    vector<Real>& COSVEC = *cptr;

    // initialize
    COSVEC[0] = 1;
    COSVEC[1] = costheta;

    //fill up COS
    for (k = 2; k<=n; k++)
    {
      if (k % 2 == 0)
      {
	COSVEC[k] = 2.0*COSVEC[k/2]*COSVEC[k/2] - 1.0;
      }
      else
      {
	COSVEC[k] = 2.0*COSVEC[k-1]*COSVEC[1]  - COSVEC[k-2];
      }
//      cout << k << "  " << COSVEC[k] << endl;
    }
    
    return;
}

struct Kprime
{
  typedef Real result_type;
  int ifrier;
  void set_ifrier(int ifr)
  {
    ifrier = ifr;
  }
  Real operator()(vec_rel_ptr_t cptr, Real t, Real u )
  {
     int k;
     Real t1, t2, t3;
     Real anum = t*t - u*u;
     Real aden = 2.0*R*(R-u);
     Real a = 2.0*R + t - u;
     Real b = 2.0*R - t - u;     
     Real den = sqrt((u+t)*a*b);

     Real Ctheta = 1.0 - anum/aden; 
     Real S = 0.0;

     cosntheta(cptr, Ctheta, ifrier); //result in cptr
     vector<Real>& C = *cptr;//(cosntheta(Ctheta, ifrier));

     t1 = 4.0*(R-u)*C[ifrier]/den;
     if (ifrier % 2 == 0) // ifrier = 0,2,4,6, ...
     {
       for (k = 1; k <= (ifrier/2); k++)
       {
	 S += C[2*k - 1];
       }
       
       t2 = -8.0*Real(ifrier)*t*t*S/(R*den);
     }
     else    // ifrier = 1,3,5,7,...
     {
       S = 0.5;
       for (k = 1; k<=((ifrier-1)/2); k++)
       {
	 S += C[2*k];
       }
       
       t2 = -8.0*t*t*Real(ifrier)*S/(R*den);
     }
 
     t3 = -2.0*t*(R-u)*(a*b - 2.0*t*(t + u))*C[ifrier]/(den*den*den);
     
     return t1+t2+t3;  

  }
};
struct GaussLegendreChangeOfVariable
{
  typedef double result_type;
  double operator()(const double x) const
  {
    return 0.25*pi*(x+1.0); 
  } 
};

template<> 
struct F_l_n<0>
{
  typedef Real result_type;
  Kprime kp;
  F_l_n<0>(int ifr)
    {
      kp.set_ifrier(ifr);
    }
  void set_ifrier(int ifr)
  {
    kp.set_ifrier(ifr);
  }
  Real operator()(vec_rel_ptr_t cptr, Real x, Real u, Real t)
  {
    double y = cos(x);    
    Real kp_arg = t*y * y  + u* ( 1 - y * y);
      if (abs(kp_arg-u) < tol)
      {
	double res  = y* y*kp(cptr, u,u);
 	return res;
      }
      else
      {

	double res  = y * y*kp(cptr, kp_arg,u);
 	return res;
      }
  }
};

template<>
struct I_l_n<0, GaussLegendreIntegration>
  {
    typedef Real result_type;
    int ifrier;
    void set_ifrier(int ifr)
    {
      ifrier = ifr;
    } 
    Real operator()(vec_rel_ptr_t cptr, Real u, Real  t)
    {
      Real scale = 2.0*sqrt(R)/(pi*sqrt(2.0*t*(R-t)));
      boost::function<Real (Real)> fa_l_n;
      fa_l_n = boost::bind(F_l_n<0>(ifrier),cptr,_1,u,t);
      #ifdef __USE_QUANTLIB__
      GaussLegendreChangeOfVariable cov;
      double integral = scale*0.25*pi*gintegrator(boost::bind(fa_l_n, boost::bind(cov, _1)));
      return integral;
      #else
      return scale*integrator(fa_l_n, 0,pi/2.0);
      #endif
    }
  };


template<int i> 
struct F_l_n
{
  typedef Real result_type;
  int ifrier;
  F_l_n<i>(int ifr)
  {
    ifrier = ifr;
  }
  Real operator()(vec_rel_ptr_t cptr, Real x, Real u, Real t)
  {
    I_l_n<i-1, AdaptiveSimponsIntegration> A;
    I_l_n<0, GaussLegendreIntegration> B;
    B.set_ifrier(ifrier);
    A.set_ifrier(ifrier);
    //    double elapsed = get_seconds();
    double l0 = A(cptr, u, x);
    double lim1 = B(cptr, x, t);
    double res = l0 * lim1;
    return res;
  }
};


template<int i, class integrator_t> 
struct I_l_n
{
typedef Real result_type;
  int ifrier;
  void set_ifrier(int ifr)
  {
    ifrier = ifr;
  }
  Real operator()(vec_rel_ptr_t cptr, Real u, Real t)
{
  boost::function<Real (Real)> fa_l_n;
  fa_l_n = boost::bind(F_l_n<i>(ifrier),cptr, _1,u,t);
  string outfn = string("L") + SSTR(i) + string("_u") + SSTR(u) + string("_t") + SSTR(t) + string(".graph");
  integrator_t integrator(0.0001, 12, outfn.c_str());
  return integrator(fa_l_n, u,t);
}
};

