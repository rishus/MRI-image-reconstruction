#include <iostream>
#define _USE_MATH_DEFINES

#include <math.h>
#include <vector>
#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
     

using std::vector;

struct itrec
{
  NDArray gquad;
  NDArray wts;
  int m_it;
  itrec(int it):m_it(it)
  {
    Matrix a_matrix = Matrix (1, 2);
    a_matrix (0) = max(4,16-floor(it/2));
    a_matrix (1) = max(4,16-floor(it/2));
    octave_value_list in = octave_value(a_matrix);
    octave_value_list out =   feval("gauss_quad_rules", in, 1);
    gquad = out(0).array_value();
    wts = out(1).array_value();
  }
};

//Nr number of iterations

struct Frec
{
  vector<double> phi;
  NDArray arg;
  NDArray wts;
  int m_nr;
  Frec(int nphi, int nr): m_nr(nr)
  {
    for(int i =0; i < nphi; ++i)
      {
	double val = (2 * M_PI * i)/nphi;
	phi.push_back(val);
      }
    for(int it = 0; it < m_nr; ++it)
      {
	
	
      }

  }
  
};


struct kernel
{
  double tol;
  int maxorder;
  double R;
  double eps;
  double Nr;
  double Nfrier;
  kernel()
  {
    tol = pow(10, -10);
    maxorder = 4;
    eps = 0.001;
    R = 1.0;
    Nr = pow(2,5);
    Nfrier = pow(2,4);
  }
  void gen_kernel()
  {
    for(int it = 0; it < Nr; ++it)
      {
	for(int nf = 0; nf  < Nfrier; ++nf)
	  {
	    Nu = all_itrecs[it].Nu;
	    for(int iu = 0; iu < Nu; ++iu)
	      {
		for(int order = 0; order < maxorder; ++order)
		  {
		    Lhandle(iu+1,order+1,ifrier+1,it+1) = add_arguments([F.u{it+1,1}(iu+1), t(it+1), ifrier, R, order]);
		  }
	      }
	  }
      }
    compute_all_integrations();
    for(int it = 0; it < Nr; ++it)
      {
	for(int nf = 0; nf  < Nfrier; ++nf)
	  {
	    Nu = all_itrecs[it].Nu;
	    for(int iu = 0; iu < Nu; ++iu)
	      {
		for(int order = 0; order < maxorder; ++order)
		  {
		    Lres(iu+1,order+1,ifrier+1,it+1) = get_res([Lhandle(iu+1,  order+1, ifrier+1,it+1)]);
		  }
	      }
	  }
      }
  }
};


int main(int  argc, char** argv)
{

  return 0;
}
