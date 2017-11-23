#include <octave/oct.h>
#include "compute_parallel.h"


DEFUN_DLD (compute_all_integrations, args, nargout,"integrate all (u,t) pairs")
{
  R = 1.0;
  compute_parallel::compute();
  double res = 0;
  return octave_value(res);
}
