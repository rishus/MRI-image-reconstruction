#include <octave/oct.h>
#include "definition.h"
#include "containers.h"
//#include "integration.h"


extern std::vector<param_t> all_arguments;
DEFUN_DLD (get_res, args, nargout,"returns result of the computation")
{
//  int nargin = args.length ();
 NDArray A = args(0).array_value();
 int idx = A(0);
 return octave_value(all_arguments[idx].res);
}

