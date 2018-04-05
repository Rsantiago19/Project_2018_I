!----------------------------------------------------------------------
! 
! 
! 
! MAIN PROGRAM FOR THE MOLECULAR DYNAMICS SIMULATION 
! 
! 
! 
!----------------------------------------------------------------------
program dynamics
! USE THE REQUIRED MODULES
use moment_module
use positions_module
use velocities_module
use andersen_therm_module
use distort_module
use vel_verlet_module
use lj_module
use pbc_module
use read_data_module
use print_positions_module
use kinetic_energy_module
use print_data_module
use send_recv_module
use rows_per_proc_module

implicit none
character(25)                           :: fName, fff
integer                                 :: i, j, k, l, m, n
integer                                 :: fStat, un
real(4)                                 :: start, finish
real(8)                                 :: dt, boxSize, cutOff, T, density
integer                                 :: nPartDim, nPart, nSteps
real(8), allocatable, dimension(:,:)    :: pos, F, vel
real(8)                                 :: V, eps, sig, time, KE, Tinst
real(8), allocatable, dimension(:)      :: total_momentum
integer                                 :: seed, trjCount, thermCount
integer                                 :: initUn, finUn, trajUn, dataUn, velUn, paramUn
integer                                 :: ierror, rank, numProcs, status, numParts
integer                                 :: myFirstPart, myLastPart
integer, parameter                      :: rMaster = 0
real(4)                                 :: printIniT,  printFinT,  printTotT
real(4)                                 :: verletIniT, verletFinT, verletTotT, verletTime
real(4)                                 :: thermoIniT, thermoFinT, thermoTotT, thermoTime
real(4)                                 :: forcesIniT, forcesFinT, forcesTotT, forcesTime
real(4)                                 :: momentIniT, momentFinT, momentTotT, momentTime
real(4)                                 :: kinetiIniT, kinetiFinT, kinetiTotT, kinetiTime


! INITIALIZE THE MPI ENVIROMENT
call mpi_init(ierror)
call mpi_comm_rank(mpi_comm_world, rank, ierror)
call mpi_comm_size(mpi_comm_world, numProcs, ierror)

! THE MASTER THREAD (rMaster) IS THE ONLY ONE OPENNING 
! AND CONTROLLING THE INPUT AND OUTPUT FILES.
if (rank == rMaster) then
        call cpu_time(start)
        call get_command_argument(1, fName, status=fStat)
        if (fStat /= 0) then
                print*, 'Any file given ---> Exitting program'
                call mpi_finalize(ierror)
                call exit()
        end if 
        un = 100
        open(unit=un, file=trim(fName), status='old') 
    
        ! READS THE DATA FROM THE INPUT FILE
        call readData(un, dt, boxSize, cutOff, nPartDim, T, eps, sig, nSteps, density, seed)
        nPart = nPartDim**3

        ! CHECK THAT THE PROVIDED NUMBER OF PARTICLES FITS THE
        ! WAY THAT THE SIMULATON BOX WILL BE BUILT
        if (mod(float(nPart),8.) /= 0.) then
                print*, 'Total number of particles must be'
                print*, 'divisible by 8 to fit a SC cell'
                print*, 'Exitting Program'
                call mpi_finalize(ierror)
                call exit()
        end if
end if

! THE MASTER THREAD BROADCAST THE DAT ATO THE REST
call mpi_bcast(dt, 1, mpi_real8, rMaster, mpi_comm_world, ierror)
call mpi_bcast(boxSize, 1, mpi_real8, rMaster, mpi_comm_world, ierror)
call mpi_bcast(cutOff, 1, mpi_real8, rMaster, mpi_comm_world, ierror)
call mpi_bcast(T, 1, mpi_real8, rMaster, mpi_comm_world, ierror)
call mpi_bcast(eps, 1, mpi_real8, rMaster, mpi_comm_world, ierror)
call mpi_bcast(sig, 1, mpi_real8, rMaster, mpi_comm_world, ierror)
call mpi_bcast(nPart, 1, mpi_integer, rMaster, mpi_comm_world, ierror)
call mpi_bcast(nSteps, 1, mpi_integer, rMaster, mpi_comm_world, ierror)
call mpi_bcast(seed, 1, mpi_integer, rMaster, mpi_comm_world, ierror)

! SEED IS MODIFIED TO GENERATE DIFFERENT RANDOM SEQUENCES AT EACH THREAD
seed = seed + rank
call srand(seed)

! SPLITS THE NUMBER OF PARTICLES THAT EVERY PROCESSOR HAVE TO WORK WITH
call rows_per_proc(nPart, myFirstPart, myLastPart)
allocate(pos(nPart,3), F(myFirstPart:myLastPart,3), vel(nPart,3), total_momentum(3))

! SETTING UP THE INITIAL CONFIGURATION AND VELOCITIES
! IT IS ONLY DONE BY THE MASTER THREAD
if (rank == rMaster) then
        ! GENERATES THE INITIAL CONFIGURATION ACCORDING TO A SIMPLE CUBIC
        ! CELL AND DISTORT IT TO BREAK THE RETICULAR ENERGY
        call SC_init_conditions(nPart, pos, boxSize)
        call distort_geometry(nPart, pos, boxSize, seed)

        ! GENERATES A VELOCITY DISTRIBUTION ACCORDING TO THE PROPER 
        ! MAXWELL-WOLTZMAN FUNCTION USING THE BOX-MULLER ALGORITHM
        call IN_velocities(nPart, T, seed, vel)
        
        ! OPEN THE OUTPUT FILES
        initUn = 101; finUn = 102; trajUn = 103; dataUn = 104; velUn = 105; paramUn = 106
        open(unit=initUn, file='initial.out')           ! Coord. Inicials
        open(unit=finUn , file='final.out')             ! Coord. Finals
        open(unit=trajUn, file='traj.xyz')              ! Trajectoria
        open(unit=dataUn, file='data.out')              ! T, Ken, V, temps...
        open(unit=velUn,  file='velocity.out')          ! T, Ken, V, temps...
        open(unit=paramUn,file='parameters.out')        ! Parameters needed for statistics
                                                ! (MB and RDF)

        write(dataUn,*) '# time, V, KE, KE + V, mom, Tins'
        write(paramUn,*) nSteps/100 - 1, nPart
        write(paramUn,*) boxSize
end if

! BROADCAST THE POSITIONS AND VELOCITIES TO ALL THE PROCESSORS
call mpi_bcast(pos, 3*nPart, mpi_real8, rMaster, mpi_comm_world, ierror)
call mpi_bcast(vel, 3*nPart, mpi_real8, rMaster, mpi_comm_world, ierror)

! FIRST PRINTS
if (rank == rMaster) then
        call print_positions(initUn, nPart, pos, time)
        call print_positions(trajUn, nPart, pos, 0.0D0)
        call print_positions(velUn,  nPart, vel, 0.0D0)
end if

! SET THE TOTAL TIMES TO ZERO 
verletTotT = 0.0
thermoTotT = 0.0
printTotT  = 0.0
forcesTotT = 0.0
momentTotT = 0.0
kinetiTotT = 0.0

! INITIALIZE THE SIMULATION TIME TO ZERO AMONG OTHER COUNTERS
time = 0.0D0; trjCount = 0; thermCount = 0

! INTEGRATE POSTIONS, CALCULATE VELOCITIES AND REESCALE WITH THE 
! THERMOSTAT AS MANY TIMES REQUIRED IN THE INPUT FILE
do i = 1, nSteps, 1

        ! INTEGRATION OF THE POSITIONS AND VELOCITIES
        call cpu_time(verletIniT)
        call vel_verlet(time, dt, pos, vel, nPart, eps, sig, boxSize, cutOff, V, F, myFirstPart, myLastPart&
        &, rank, status, forcesTotT)
        call cpu_time(verletFinT)
        verletTotT = verletTotT + (verletFinT - verletIniT)

        ! CALCULATE THE MOMENTUM AND THE KINETIC ENERGY AND PRINT DATA
        if (mod(trjCount,100) == 0) then
                call cpu_time(momentIniT)
                call momentum(myFirstPart, myLastPart, vel, total_momentum)
                call cpu_time(momentFinT)
                momentTotT = momentTotT + (momentFinT - momentIniT)

                call cpu_time(kinetiIniT)
                call kinetic_energy(vel, KE, Tinst, myFirstPart, myLastPart, nPart, rank)
                call cpu_time(kinetiFinT)
                kinetiTotT = kinetiTotT + (kinetiFinT - kinetiIniT)

                ! PRINT THE DATA
                if (rank == rMaster) then
                        call cpu_time(printIniT)
                        call print_positions(trajUn, nPart, pos, time)
                        call print_positions(velUn,  nPart, vel, time)
                        call print_data(time, V, KE, Tinst, total_momentum, dataUn)
                        call cpu_time(printFinT)
                        printTotT = printTotT + (printFinT - printIniT)
                end if
                trjCount = 0
        end if

        ! REESCALE THE VELOCITIES USING THE ANDERSEN THERMOSTAT
        if (mod(thermCount,10) == 0) then
                call cpu_time(thermoIniT)
                call andersen_thermo(dt, T, nPart, seed, vel, myFirstPart, myLastPart, rank, status)
                call cpu_time(thermoFinT)
                thermoTotT = thermoTotT + (thermoFinT - thermoIniT)
                thermCount = 0
        end if
        trjCount = trjCount + 1
        thermCount = thermCount + 1
end do

! TAKE THE TOTAL CPU-TIME'S 
call mpi_reduce(verletTotT, verletTime, 1, mpi_real4, mpi_sum, rMaster, mpi_comm_world, ierror) 
call mpi_reduce(thermoTotT, thermoTime, 1, mpi_real4, mpi_sum, rMaster, mpi_comm_world, ierror) 
call mpi_reduce(forcesTotT, forcesTime, 1, mpi_real4, mpi_sum, rMaster, mpi_comm_world, ierror) 
call mpi_reduce(momentTotT, momentTime, 1, mpi_real4, mpi_sum, rMaster, mpi_comm_world, ierror) 
call mpi_reduce(kinetiTotT, kinetiTime, 1, mpi_real4, mpi_sum, rMaster, mpi_comm_world, ierror) 

! TREAT, MODIFY AND PRINT CPU-TIMES TO PRINT THE CPU-TIME PER PROCESSOR
if (rank == rMaster) then
        call print_positions(finUn, nPart, pos, time)
        call cpu_time(finish)
        verletTime = verletTime - forcesTime
        write(dataUn,*) '# VERLET TOTAL_CPU_TIME  : ', verletTime, 'CPU_TIME_PER_PROC: ', verletTime/numProcs
        write(dataUn,*) '# THERMO TOTAL_CPU_TIME  : ', thermoTime, 'CPU_TIME_PER_PROC: ', thermoTime/numProcs
        write(dataUn,*) '# FORCES TOTAL_CPU_TIME  : ', forcesTime, 'CPU_TIME_PER_PROC: ', forcesTime/numProcs
        write(dataUn,*) '# MOMENT TOTAL_CPU_TIME  : ', momentTime, 'CPU_TIME_PER_PROC: ', momentTime/numProcs
        write(dataUn,*) '# KINETIC TOTAL_CPU_TIME : ', kinetiTime, 'CPU_TIME_PER_PROC: ', kinetiTime/numProcs
        write(dataUn,*) '# PRINTING TIME (MASTER) : ', printTotT
        write(dataUn,*) "# TOTAL CPU TIME         : ", finish - start
        write(dataUn,*)
        write(dataUn,*)
        write(dataUn,*) '# RELATIVE TIME TO THE TOTAL CPU_TIME (IN %):'
        write(dataUn,*) '# VERLET   : ', 100*verletTime/(numProcs*(finish-start))
        write(dataUn,*) '# THERMO   : ', 100*thermoTime/(numProcs*(finish-start))
        write(dataUn,*) '# FORCES   : ', 100*forcesTime/(numProcs*(finish-start))
        write(dataUn,*) '# MOMENT   : ', 100*momentTime/(numProcs*(finish-start))
        write(dataUn,*) '# KINETIC  : ', 100*kinetiTime/(numProcs*(finish-start))
        write(dataUn,*) '# PRINTING : ', 100*printTotT/(numProcs*(finish-start))

        close(un); close(initUn); close(finUn); close(trajUn); close(dataUn)
end if

call mpi_finalize(ierror)
contains

end program dynamics
    
