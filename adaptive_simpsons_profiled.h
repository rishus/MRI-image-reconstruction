

#include <math.h>
#include <boost/math/distributions.hpp>
#include <boost/function.hpp>
#include <fstream>

#ifndef __ADAPTIVESIMPSONINTEGRATIONPROFILED_H
#define __ADAPTIVESIMPSONINTEGRATIONPROFILED_H
class AdaptiveSimponsIntegrationProfiled
{

public:
  std::ofstream ofs;
Real m_epsilon;// = 0.000000001;
int m_maxRecursionDepth;// = 15;
 AdaptiveSimponsIntegrationProfiled(Real eps = 0.0001, int dep=12, const char* fn  = (char*) 0):m_epsilon(eps),m_maxRecursionDepth(dep)
{
  ofs.open(fn);

}

Real adaptiveSimpsonsAux(boost::function<Real (Real)>& integrand, Real a, Real b, Real epsilon,                 
                         Real S, Real fa, Real fb, Real fc, int bottom) 
{                 
Real c = (a + b)/2, h = b - a;
Real d = (a + c)/2, e = (c + b)/2;                                                              
Real fd = integrand(d), fe = integrand(e);                                                                      
Real Sleft = (h/12)*(fa + 4*fd + fc);                                                           
Real Sright = (h/12)*(fc + 4*fe + fb);                                                          
Real S2 = Sleft + Sright;                                                                       
if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon) 
  {
    //std::cout<<c-d<<" "<<"["<<c<<" "<<d<<"]="<<S<<std::endl;
    //     //std::cout<<"adaptivesimpsons depth"<<m_maxRecursionDepth - bottom<<" accuracy = "<<fabs(S2-S)<<std::endl;
    ofs<<a<<" "<<b<<" "<<S<<std::endl;
    return S2 + (S2 - S)/15;                                                                        
  }

return adaptiveSimpsonsAux(integrand, a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1) +                    
adaptiveSimpsonsAux(integrand, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1);                     
}         
 

Real adaptiveSimpsons(boost::function<Real (Real)>& integrand,   // ptr to function
                           Real a, Real b,  // interval [a,b]
                           Real epsilon,  // error tolerance
                           int maxRecursionDepth) {   // recursion cap        
  Real c = (a + b)/2, h = b - a;                                                                  
  Real fa = integrand(a), fb = integrand(b), fc = integrand(c);                                                           
  Real S = (h/6)*(fa + 4*fc + fb);                                                                
return adaptiveSimpsonsAux(integrand, a, b, epsilon, S, fa, fb, fc, maxRecursionDepth);                   
}        

Real operator()(boost::function<Real (Real)>& integrand, Real a, Real b)
{
return adaptiveSimpsons(integrand, a, b, m_epsilon, m_maxRecursionDepth);
}


};

#endif
