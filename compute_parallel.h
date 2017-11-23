#include <stdlib.h>
#include <algorithm>
#include "definition.h"
#include "integration.h"
#include "containers.h"
#include <string>

#include <fstream>
#include <sstream>
#include <algorithm>
#include "adaptive_simpsons.h"
#include "adaptive_simpsons_profiled.h"


#include "mpi.h"
extern std::vector<param_t> all_arguments;
std::vector<int> work_queue;

namespace compute_parallel
{
  #define num_threads 12
  void* thread_entry(void* arg)
  {
    int i ;
    long tid;
    tid = (long)arg;
    int num_args = all_arguments.size();
    int range = num_args/num_threads;
    int lidx = range * tid; 
    string sstid = SSTR(tid);

    string fn = string("progress_") + sstid  + string(".dat");
    std::ofstream ofs(fn.c_str());
	 
    int uidx =  range * (tid+1);
    std::cerr<<"range "<<tid<<" "<<lidx<<" "<<uidx<<std::endl;
    I_l_n<0, GaussLegendreIntegration> L0;
    I_l_n<1, AdaptiveSimponsIntegrationProfiled> L1;
    I_l_n<2, AdaptiveSimponsIntegrationProfiled> L2;
    I_l_n<3, AdaptiveSimponsIntegrationProfiled> L3;
    I_l_n<4, AdaptiveSimponsIntegrationProfiled> L4;
    vec_rel_ptr_t cptr=  new vector<Real>(1000,0.0);
    if(uidx > num_args) uidx = num_args;
    double elapsed;
    for(int  i = lidx;  i < uidx; ++i)
      {
	param_t& arg = all_arguments[work_queue[i]];
	Real u = arg.u;
	Real t= arg.t;
	int ifrier = arg.ifrier;
	int l = arg.level;
	double res;
	elapsed = get_seconds();
	//call appropriate Li
	if(l == 0){ L0.set_ifrier(ifrier); res = L0(cptr, u,t);}
	if(l == 1){ L1.set_ifrier(ifrier); res = L1(cptr, u,t);}
	if(l == 2){ L2.set_ifrier(ifrier); res = L2(cptr, u,t);}
	if(l == 3){ L3.set_ifrier(ifrier); res = L3(cptr, u,t);}
	if(l == 4){ L4.set_ifrier(ifrier); res = L4(cptr, u,t);}
	arg.res = res;
	string outfn = string("L") + SSTR(i) + string("_u") + SSTR(u) + string("_t") + SSTR(t) + string(".graph");
	ofs<<"handle="<<i<<" u="<<u<<" t="<<t<<" ifrier="<<ifrier<<" level="<<l<<" time="<<get_seconds() - elapsed<<" "<<outfn<<endl;
	if((i - lidx)%50 == 0)
		fprintf(stderr, "tid = %d, %ld/%ld\n", tid, (i-lidx), (uidx-lidx));
      }
    delete cptr;
    pthread_exit((void*) arg);
  }

  void compute()
  {
    MPI_Init () 
    //run stuff and wait till doneo
    pthread_t  all_threads[num_threads];
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    void *status;      
    assert(all_arguments.size() >= num_threads);
    for(int i = 0; i < all_arguments.size(); ++i)
	{
		work_queue.push_back(i);
	}
    std::random_shuffle(work_queue.begin(), work_queue.end());
    for(long i =0; i < num_threads; ++i)
      {
	pthread_create(&all_threads[i], &attr, thread_entry, (void*)i);
      }
    pthread_attr_destroy(&attr); 
    for(int t=0; t < num_threads; t++) {
      pthread_join(all_threads[t], &status);
    }
  }


}



