module andersen_therm_module
use mpi
use send_recv_module
contains
subroutine andersen_thermo(dt, T, nPart, seed, vel, myFirstPart, myLastPart, rank, status)
implicit none
integer, intent(in)                             :: nPart, seed, myFirstPart, myLastPart
integer, intent(in)                             :: rank, status
real(8), intent(in)                             :: dt, T
real(8), dimension(nPart,3), intent(inout)      :: vel
real(8)                                         :: x1, x2, y1, y2, w
integer                                         :: i, j, k

! EACH THREAD ONLY OPERATES TO ITS PARTICLE RANGE
! SO ONLY LOOPS FROM myFirstPart TO myLastPart
do i = myFirstPart, myLastPart, 1
        if (rand() < 0.1) then
                ! canvia la velocitat d'aquesta particula segons una distribució
                ! gaussiana com maxwell boltzman.
                ! Aixó ho fem utilitzant el algoritme Box_muller, pero com
                ! aquest genera parells de numeros aleatoris i tenim tres
                ! components en el vector de velocitats, podem fer dues vegades
                ! el BM i descartar un dels numeros aleatoris.
                w = 1.0D0
                do while (w >= 1.0D0)
                        x1 = 2.0*rand() - 1.0
                        x2 = 2.0*rand() - 1.0
                        w = x1**2. + x2**2.
                end do
                w  = dsqrt(-2.0*dlog(w)/w)
                y1 = x1*w
                y2 = x2*w
                vel(i,1) = y1*dsqrt(T)
                vel(i,2) = y2*dsqrt(T)
                w = 1.0D0
                do while (w >= 1.0D0)
                        x1 = 2.0*rand() - 1.0
                        x2 = 2.0*rand() - 1.0
                        w = x1**2. + x2**2.
                end do
                w  = dsqrt(-2.0*dlog(w)/w)
                y1 = x1*w
                vel(i,3) = y1*dsqrt(T)
        end if
end do
call send_recv_array(vel(myFirstPart:myLastPart,:),myFirstPart,myLastPart,rank,nPart,status,vel)
end subroutine andersen_thermo
end module andersen_therm_module
