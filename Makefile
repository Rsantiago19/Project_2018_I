compile-pbc: pbc_module.f90
	gfortran -c pbc_module.f90

compile-lj: lj_module.f90
	gfortran -c lj_module.o lj_module.f90


