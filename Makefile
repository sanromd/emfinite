
fdtd_da_pml.so: fdtd_da_pml.f90
	f2py --f90flags='-g3 -fcheck=all' -m $(basename $(notdir $@)) -c $^
