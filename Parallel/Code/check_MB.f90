program check_MB
use mpi
use rows_per_proc_module
implicit none
character(25)                           :: fName
character(5)                            :: trash
integer                                 :: fStat, un, unOut, paramUn
integer                                 :: nPart, nIt, nHist
integer                                 :: i, j, k, l
real(8), allocatable, dimension(:,:)    :: vel
real(8), allocatable, dimension(:)      :: histogram, masterHistogram
real(8), dimension(3)                   :: vec
real(8)                                 :: iniHist, finHist, pasH, modV
real(8)                                 :: minVel, minVel2
real(8)                                 :: iniT, finalT, temps
!Variables MPI
integer                                 :: ierror, rank, numProcs, status, numParts, myFirstPart, myLastPart
integer, parameter                      :: rMaster = 0

! INITIALIZE ALL THE MPI ENVIROMENT
call mpi_init(ierror)
call mpi_comm_rank(mpi_comm_world, rank, ierror)
call mpi_comm_size(mpi_comm_world, numProcs, ierror)


! ALL THE THREADS READ THE INPUT THAT SHOULD BE PROVIDED
! OTHERWISE THE PROGRAM ENDS BY CALLING mpi_finalize()
call get_command_argument(1, fName, status=fStat)
if (fStat /= 0) then
        print*, 'Any file given ---> Exitting program'
        call mpi_finalize(ierror)
        call exit()
end if
un = 100; unOut = 101
open(unit=un, file=trim(fName), status='old')
call cpu_time(iniT)
call get_command_argument(2, fName, status=fStat)
if (fStat /= 0) then
        print*, 'Any file given ---> Exitting program'
        call mpi_finalize(ierror)
        call exit()
end if
open(unit=paramUn, file=trim(fName), status='old')
read(paramUn,*) nIt, nPart


! SET THE MINIMUM AND MAXIMUM VELOCITIES CATCHED TO -15:15 
! AND CALCULATE THE HISTOGRAM WIDTH
iniHist = -15.0D0; finHist = 15.0D0; nHist = 3000
pasH = (finHist - iniHist)/dfloat(nHist)
call mpi_bcast(nIt, 1, mpi_integer, rMaster, mpi_comm_world, ierror)
call mpi_bcast(nPart, 1, mpi_integer, rMaster, mpi_comm_world, ierror)
allocate(vel(nPart,3), histogram(nHist + 2), masterHistogram(nHist + 2))

! SPLITS THE PARTICLES THAT EVERY PROCESSOR WILL WORK WITH AMONG THE THREADS
call rows_per_proc(nPart, myFirstPart, myLastPart)

! EVERY THREAD FILLS THE PARTIAL HISTOGRAM
histogram(:) = 0
masterHistogram(:) = 0
do i = 1, nIt, 1
        print *, i
        read(un,*) nPart
        read(un,*) trash
        do j = 1, nPart, 1
                read(un,*) trash, vel(j,:)
        end do
        do k = myFirstPart, myLastPart, 1; do l = 1, 3, 1
                vec(:) = vel(k,:)

                minVel  = iniHist
                minVel2 = iniHist + pasH
                do j = 1, nHist, 1
                        if ((vec(l) <= minVel2).and.(vec(l) > minVel)) then
                                histogram(j+1) = histogram(j+1) + 1
                        end if
                        minVel  = minVel  + pasH
                        minvel2 = minVel2 + pasH
                end do

                if (vec(l) < iniHist) then
                        histogram(1) = histogram(1) + 1
                else if (modV >= finHist) then
                        histogram(nHist+2) = histogram(nHist+2) + 1
                end if
        end do; end do
end do

! THE TOTAL HISTOGRAM IS COLLECTED BY THE MASTER
call mpi_reduce(histogram, masterHistogram, (nHist + 2), mpi_real8, mpi_sum, rMaster, mpi_comm_world, ierror)

if (rank == rMaster) then
        ! NORMALIZATION OF THE HISTOGRAM AND OUTPUT O THE FILE
        masterHistogram(:) = masterHistogram(:)/sum(3*masterHistogram(:)*pasH)

        open(unit=unOut, file='MB_velocityDistribution.out')
        do i = 1, nHist + 2, 1
                write(unOut,*) iniHist + (i-1)*pasH, masterHistogram(i)
        end do
        close(un); close(unOut); close(paramUn)
        call cpu_time(finalT)
        temps = finalT - iniT
        print *, "CPU_TIME:", temps
end if

call mpi_finalize(ierror)

end program check_MB
