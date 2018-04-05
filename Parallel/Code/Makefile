SHELL:=/bin/bash
path = $(shell pwd)/../
dir_name= test_run
new_dir = $(path)$(dir_name)
num_procs = 2
all_exec = compile-main compile-MB compile-rdf
to_compile = compile-rows_per_proc compile-send_array compile-pbc compile-lj compile-mom compile-init_pos compile-init_vel compile-thermo compile-verlet compile-pos compile-KE compile-read_data compile-print_data compile-distort 
.PHONY: clean backup run-main run-MB run-rdf variables module-load

compile-all: $(all_exec) 

module-load:
	source $(MODULEPATH)/init/bash \
	module load openmpi

compile-main: EIA_dynamics_main.f90 $(to_compile)
	mpif90 -g -o 	dynamics_mpi \
			send_recv_module.o \
			pbc_module.o \
			distort_module.o \
		       	lj_module.o \
			moment_module.o \
			position_module.o \
			rows_per_proc_module.o \
			velocities_module.o \
			andersen_therm_module.o \
			vel_verlet_module.o \
			print_positions_module.o \
			kinetic_energy_module.o \
			read_data_module.o \
			print_data_module.o \
			EIA_dynamics_main.f90
			
compile-MB: compile-rows_per_proc check_MB.f90
	mpif90 -g -o 	check_MB_mpi \
			rows_per_proc_module.o \
			check_MB.f90

compile-rdf: compile-rows_per_proc rdf.f90 
	mpif90 -g -o 	rdf_mpi \
			rows_per_proc_module.o \
			rdf.f90 \
			pbc_module.f90

compile-pbc: pbc_module.f90 
	mpif90 -c -g	 pbc_module.f90

compile-lj: lj_module.f90
	mpif90 -c -g 	lj_module.f90

compile-mom: moment_module.f90
	mpif90 -c -g	moment_module.f90

compile-init_pos: position_module.f90
	mpif90 -c -g 	position_module.f90

compile-init_vel: velocities_module.f90
	mpif90 -c -g  	velocities_module.f90

compile-thermo: andersen_therm_module.f90
	mpif90 -c -g  	andersen_therm_module.f90

compile-verlet: vel_verlet_module.f90
	mpif90 -c -g  	vel_verlet_module.f90

compile-pos: print_positions_module.f90
	mpif90 -c -g  	print_positions_module.f90

compile-KE: kinetic_energy_module.f90
	mpif90 -c -g  	kinetic_energy_module.f90

compile-read_data: read_data_module.f90
	mpif90 -c -g  	read_data_module.f90

compile-print_data: print_data_module.f90
	mpif90 -c -g  	print_data_module.f90

compile-distort: distort_module.f90
	mpif90 -c -g  	distort_module.f90

compile-send_array: send_recv_module.f90
	mpif90 -c -g  	send_recv_module.f90

compile-rows_per_proc: rows_per_proc_module.f90
	mpif90 -c -g 	rows_per_proc_module.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f dynamics_mpi rdf_mpi check_MB_mpi
	rm -f *.out traj.xyz

backup:
	@if [ -d 'backup' ]; then \
		rm -f ./backup/*; \
	else \
		mkdir ./backup; \
	fi 
	cp *.f90 ./backup

run-main: compile-main
	@if [ ! -d $(new_dir) ]; then \
		mkdir $(new_dir); \
		cp ./input.dat $(new_dir); \
		cp ./dynamics_mpi $(new_dir); \
		cd $(new_dir); \
		mpirun -np $(num_procs) dynamics_mpi input.dat; \
	else \
		echo 'THE DIRECTORY ALREADY EXISTS'; \
		echo 'PLEASE CHOOSE AN OTHER NAME' ; \
	fi 

run-main-bsc: compile-main
	@if [ ! -d $(new_dir) ]; then \
                mkdir $(new_dir); \
                cp ./input.dat $(new_dir); \
                cp ./dynamics_mpi $(new_dir); \
		cp ./sub_bsc.sh $(new_dir); \
                cd $(new_dir); \
                bsub -W 210 -q training sub_bsc.sh $(num_procs); \
        else \
                echo 'THE DIRECTORY ALREADY EXISTS'; \
                echo 'PLEASE CHOOSE AN OTHER NAME' ; \
        fi


run-MB: compile-MB
	@if [ -d $(new_dir) ]; then \
		cp ./check_MB_mpi $(new_dir); \
		cd $(new_dir); \
		mpirun -np $(num_procs) check_MB_mpi velocity.out parameters.out; \
	else \
		echo 'THE DIRECTORY DOES NOT EXIST'; \
		echo 'PLEASE TARGET A VALID ONE'; \
	fi

run-rdf: compile-rdf
	@if [ -d $(new_dir) ]; then \
		cp ./rdf_mpi $(new_dir); \
		cd $(new_dir); \
		mpirun -np $(num_procs) rdf_mpi traj.xyz parameters.out; \
	else \
		echo 'THE DIRECTORY DOES NOT EXIST'; \
		echo 'PLEASE TARGET A VALID ONE'; \
	fi

help: 
	@echo '-------------------------------------------------'
	@echo Type "make" to compile everything
	@echo Type "make run-main" to run the MD simulation in ../test_dir and 2 threads by default.
	@echo '        'To change default values: "make run-main num_procs=# dir_name='your/dir/'"
	@echo Type "make run-rdf" to run the RDF computation program.
	@echo Type "make run-MB" to perform the checking the velocity distribution.
	@echo '        '\(This two commands must be executed AFTER running MD.
	@echo '        'But the Makefile already have this feature into account.
	@echo '        'Default options for "run-main" apply to "make run-rdf" and "make run-MB" too.\)
	@echo Type "make clean" to remove *.o, *.mod, *.xyz, *.out and the exec files.
	@echo Type "make backup" to create a backup directory with the *.f90 files.
	@echo Type "make help" to show the aviable options.
	@echo 
	@echo 
	@echo 
	@echo Default options to the "run" commands. 
	@echo 
	@echo 
	@echo "'dir_name' default is "$(dir_name)
	@echo 
	@echo 
	@echo "'num_procs' default is "$(num_procs)
	@echo 
	@echo '-------------------------------------------------'




























