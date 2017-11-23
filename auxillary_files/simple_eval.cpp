#include <iostream>
#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
     
int
main (void)
{
  string_vector argv (2);
  argv(0) = "embedded";
  argv(1) = "-q";
     
  octave_main (2, argv.c_str_vec(), 1);
     
  octave_idx_type n = 2;
  Matrix a_matrix = Matrix (1, 2);
  a_matrix (0) = 8;
  a_matrix (1) = 8;
  

  octave_value_list in = octave_value(a_matrix);
  octave_value_list out =   feval("gauss_quad_rules", in, 1);
  NDArray x = out(0).array_value();
  NDArray y = out(1).array_value();

  std::cout<< x.dims()(0)<<"," <<x.dims()(1)<<std::endl;
  std::cout<< y.dims()(0)<<"," <<y.dims()(1)<<std::endl;
  std::cout<< x <<std::endl;
  std::cout<< y <<std::endl;
  return 0;
}
