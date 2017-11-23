all: containers.so add_arguments.oct get_res.oct compute_all_integrations.oct

containers.so : containers.cpp containers.h
	g++ -fPIC -shared containers.cpp -o libcontainers.so


add_arguments.oct : add_arguments.cpp 
	mkoctfile add_arguments.cpp  -lQuantLib -lcontainers -L.

get_res.oct : get_res.cpp 
	mkoctfile get_res.cpp -lQuantLib -lcontainers -L.

clear_arguments.oct : clear_arguments.cpp 
	mkoctfile clear_arguments.cpp -lQuantLib -lcontainers -L.

compute_all_integrations.oct : compute_all_integrations.cpp compute_parallel.h integration.h
	mkoctfile compute_all_integrations.cpp -lQuantLib -lcontainers -L.

clean:
	rm add_arguments.oct containers.so get_res.oct compute_all_integrations.oct
