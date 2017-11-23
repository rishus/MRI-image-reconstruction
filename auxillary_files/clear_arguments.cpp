#include <octave/oct.h>
#include "definition.h"
#include "containers.h"
//#include "integration.h"


extern std::vector<param_t> all_arguments;
DEFUN_DLD (clear_arguments, args, nargout,"clear the vector of arguments")
{
//  int nargin = args.length ();
  all_arguments.clear();
  double res = 0.0;
  return octave_value(res);
}

