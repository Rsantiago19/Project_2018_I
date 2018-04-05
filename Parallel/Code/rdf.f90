program rdf
use pbc_module
use mpi
use rows_per_proc_module
implicit none
character(25)                           :: fName
character(5)                            :: trash
real(8)                                 :: boxSize
integer                                 :: fStat, un, unOut, paramUn
integer                                 :: nPart, nIt, nRad
integer                                 :: i, j, k, l
real(8), allocatable, dimension(:,:)    :: pos
real(8), allocatable, dimension(:)      :: histogram, masterHistogram
real(8), dimension(3)                   :: vec, tar
real(8)                                 :: iniRad, finRad, pasR, modV
real(8)                                 :: minRad, minRad2, factor
real(8)                                 :: iniT, finalT, temps, Vmin, Vmax
!Variables MPI
integer                                 :: ierror, rank, numProcs, status, numParts, myFirstPart, myLastPart
integer, parameter                      :: rMaster = 0

call mpi_init(ierror)
call mpi_comm_rank(mpi_comm_world, rank, ierror)
call mpi_comm_size(mpi_comm_world, numProcs, ierror)
call get_command_argument(1, fName, status=fStat)
if (fStat /= 0) then
        print*, 'Any file given ---> Exitting program'
        call mpi_finalize(ierror)
        call exit()
end if

un = 100; unOut = 101
open(unit=un, file=trim(fName), status='old')

if (rank == rMaster) then
        call cpu_time(iniT)
        call get_command_argument(2, fName, status= fStat)
        if (fStat /= 0) then
                print*, 'Size of the box needed ---> Exitting program'
                call mpi_finalize(ierror)
                call exit()
        end if
        open(unit=paramUn, file=trim(fName), status='old')
        read(paramUn,*) nIt, nPart
        read(paramUn,*) boxSize

end if

call mpi_bcast(nPart, 1, mpi_integer, rMaster, mpi_comm_world, ierror)
call mpi_bcast(nIt, 1, mpi_integer, rMaster, mpi_comm_world, ierror)
call mpi_bcast(boxSize, 1, mpi_real8, rMaster, mpi_comm_world, ierror)

iniRad = 1D-4; finRad = boxSize; nRad = 200
pasR = (finRad - iniRad)/dfloat(nRad)

allocate(pos(nPart,3), histogram(nRad + 2), masterHistogram(nRad + 2))

call rows_per_proc(nPart, myFirstPart, myLastPart)

histogram(:) = 0.0D0
masterHistogram(:) = 0.0D0
do i = 1, nIt, 1
        if (rank == rMaster) print*, 'jejej he fet una iter.', i
        read(un,*) nPart
        read(un,*) trash
        do j = 1, nPart, 1
                read(un,*) trash, pos(j,:)
        end do
 
        do l = myFirstPart, myLastPart, 1
                tar(:) = pos(l,:)
                do k = 1, nPart, 1
                        if (k == l) cycle
                        vec(:) = pos(k,:) - tar(:)
                        call pbc(vec, boxSize)
                        modV = dsqrt(dot_product(vec, vec))

                        minRad  = iniRad
                        minRad2 = iniRad + pasR
                        do j = 1, nRad, 1
                                if ((modV <= minRad2).and.(modV > minRad)) then
                                        histogram(j+1) = histogram(j+1) + 1
                                end if
                                minRad  = minRad2
                                minRad2 = minRad2 + pasR
                        end do

                        if (modV < iniRad) then
                                histogram(1) = histogram(1) + 1
                        else if (modV >= finRad) then
                                histogram(nRad+2) = histogram(nRad+2) + 1
                        end if
                end do
        end do
end do

call mpi_reduce(histogram, masterHistogram, (nRad + 2), mpi_real8, mpi_sum, rMaster, mpi_comm_world, ierror)


!histogram(:) = histogram(:)/(pasR*sum(histogram))
if (rank == rMaster) then
        open(unit=unOut, file='rdf.out')
        minRad  = iniRad
        minRad2 = iniRad + pasR
        do i = 1, nRad + 2, 1
                Vmin = 4./3.*3.14*minRad**3
                Vmax = 4./3.*3.14*minRad2**3
                factor = Vmax - Vmin
                write(unOut,*) iniRad + (i-1)*pasR, masterHistogram(i)/(nIt*factor*nPart)
                minRad  = minRad2
                minRad2 = minRad2 + pasR
        end do
        close(un); close(unOut); close(paramUn)
        
        call CPU_TIME(finalT)
        temps = finalT - iniT
        print *, "CPU_TIME:", temps
end if

call mpi_finalize(ierror)

end program rdf

