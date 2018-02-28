program rdf
use pbc_module
implicit none
character(25)                           :: fName
character(5)                            :: trash
real(8)                                 :: boxSize
integer                                 :: fStat, un, unOut
integer                                 :: nPart, nIt, nRad
integer                                 :: i, j, k, l
real(8), allocatable, dimension(:,:)    :: pos
real(8), allocatable, dimension(:)      :: histogram
real(8), dimension(3)                   :: vec, tar
real(8)                                 :: iniRad, finRad, pasR, modV
real(8)                                 :: minRad, minRad2, factor
integer                                 :: targetP

call get_command_argument(1, fName, status=fStat)
if (fStat /= 0) then
        print*, 'Any file given ---> Exitting program'
        call exit()
end if
un = 100; unOut = 101
open(unit=un, file=trim(fName), status='old')
call get_command_argument(2, fName, status= fStat)
if (fStat /= 0) then
        print*, 'Size of the box needed ---> Exitting program'
        call exit()
end if
read(fName,*) boxSize

iniRad = 0.0D0; finRad = boxSize; nRad = 200
pasR = (finRad - iniRad)/dfloat(nRad)

read(un,*) nIt, nPart
allocate(pos(nPart,3), histogram(nRad + 2))

histogram(:) = 0
do i = 1, nIt, 1
        read(un,*) nPart
        read(un,*) trash
        do j = 1, nPart, 1
                read(un,*) trash, pos(j,:)
        end do
        targetP = 1
        do l = 1, nPart, 1
        tar(:) = pos(targetP,:)
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
                                minRad  = minRad  + pasR
                                minRad2 = minRad2 + pasR
                        end do

                        if (modV < iniRad) then
                                histogram(1) = histogram(1) + 1
                        else if (modV >= finRad) then
                                histogram(nRad+2) = histogram(nRad+2) + 1
                        end if
                end do
                targetP = targetP + 1
        end do
end do

!histogram(:) = histogram(:)/(pasR*sum(histogram))

open(unit=unOut, file='rdf.out')
minRad  = iniRad
minRad2 = iniRad + pasR
factor = 4*3.14*pasR*nPart
do i = 1, nRad + 2, 1
        write(unOut,*) iniRad + (i-1)*pasR, histogram(i)/(nIt*factor*minRad**2)
        minRad = minRad + pasR
end do
close(un); close(unOut)
contains

end program rdf

