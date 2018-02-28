program check_MB
implicit none
character(25)                           :: fName
character(5)                            :: trash
integer                                 :: fStat, un, unOut
integer                                 :: nPart, nIt, nHist
integer                                 :: i, j, k, l
real(8), allocatable, dimension(:,:)    :: vel
real(8), allocatable, dimension(:)      :: histogram
real(8), dimension(3)                   :: vec
real(8)                                 :: iniHist, finHist, pasH, modV
real(8)                                 :: minVel, minVel2

call get_command_argument(1, fName, status=fStat)
if (fStat /= 0) then
        print*, 'Any file given ---> Exitting program'
        call exit()
end if
un = 100; unOut = 101
iniHist = -15.0D0; finHist = 15.0D0; nHist = 3000
pasH = (finHist - iniHist)/dfloat(nHist)

open(unit=un, file=trim(fName), status='old')
read(un,*) nIt, nPart
allocate(vel(nPart,3), histogram(nHist + 2))


histogram(:) = 0
do i = 1, nIt, 1
        read(un,*) nPart
        read(un,*) trash
        do j = 1, nPart, 1
                read(un,*) trash, vel(j,:)
        end do
        do k = 1, nPart, 1; do l = 1, 3, 1
                vec(:) = vel(k,:)
                !modV = dsqrt(dot_product(vec,vec))

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
! NORMALITZACIÓ
histogram(:) = histogram(:)/sum(3*histogram(:)*pasH)

open(unit=unOut, file='MB_velocityDistribution.out')
do i = 1, nHist + 2, 1
        write(unOut,*) iniHist + (i-1)*pasH, histogram(i)
end do
close(un); close(unOut)


contains

end program check_MB
