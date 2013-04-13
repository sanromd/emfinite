
fdtd_da_pml.so: fdtd_da_pml.f90
	f2py  -m $(basename $(notdir $@)) -c $^
