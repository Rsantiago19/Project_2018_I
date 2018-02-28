compile-pbc: pbc_module.f90
	gfortran -c pbc_module.f90

compile-lj: lj_module.f90
	gfortran -c lj_module.f90

compile-mom: moment_module.f90
	gfortran -c moment_module.f90

compile-init_pos: position_module.f90
	gfortran -c position_module.f90

compile-init_vel: velocities_module.f90
	gfortran -c velocities_module.f90

compile-thermo: andersen_therm_module.f90
	gfortran -c andersen_therm_module.f90

compile-verlet: vel_verlet_module.f90
	gfortran -c vel_verlet_module.f90

compile-pos: print_positions_module.f90
	gfortran -c print_positions_module.f90

compile-KE: kinetic_energy_module.f90
	gfortran -c kinetic_energy_module.f90

compile-read_data: read_data_module.f90
	gfortran -c read_data_module.f90

compile-print_data: print_data_module.f90
	gfortran -c print_data_module.f90

compile-main: EIA_dynamics_main.f90
	gfortran pbc_module.o lj_module.o moment_module.o position_module.o velocities_module.o andersen_therm_module.o vel_verlet_module.o print_positions_module.o kinetic_energy_module.o read_data_module.o print_data_module.o EIA_dynamics_main.f90 -o dynamics

compile-MB: check_MB.f90
	gfortran check_MB.f90 -o check_MB
compile-rdf: rdf.f90
	gfortran pbc_module.o rdf.f90 -o rdf

