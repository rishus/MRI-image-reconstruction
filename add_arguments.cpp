#include <octave/oct.h>
#include "definition.h"
#include "containers.h"
//#include "integration.h"


extern std::vector<param_t> all_arguments;
DEFUN_DLD (add_arguments, args, nargout,"add L0(u,t) work")
{
//  int nargin = args.length ();
 NDArray A = args(0).array_value();
 
 Real u = A(0);
 Real t= A(1); 
 int ifrier = A(2);
 double  R=A(3);
 Real ll = A(4);

 Real res=0.0;
 param_t p(ifrier, u, t, ll);
 
 all_arguments.push_back(p);
 int idx = all_arguments.size() - 1;
 return octave_value(idx);
}

