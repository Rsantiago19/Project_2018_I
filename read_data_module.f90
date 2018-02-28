module read_data_module
implicit none
contains

subroutine readData(un, dt, boxSize, cutOff, nPartDim, T, eps, sig, nSteps, density, seed)
implicit none
real(8), intent(out)                            :: dt, boxSize, cutOff, T
real(8), intent(out)                            :: eps, sig
integer, intent(out)                            :: nPartDim, un, nSteps, seed
real(8)                                         :: density
read(un,*) dt                   ! Increment de temps
read(un,*) density              ! Densitat en unitats adimensionals 
read(un,*) cutOff               ! Distància màxima a la que es calcula LJ
read(un,*) nPartDim             ! Numero total de partícules en cada dimensio
read(un,*) T                    ! Temperatura de la simulació
read(un,*) eps                  ! Epsilon per al potencial LJ
read(un,*) sig                  ! Sigma per al potencial LJ
read(un,*) nSteps               ! Numero de passos de temps que fa la simulació
read(un,*) seed                 ! llavor per la generació de numeros aleatoris.
close(un)
! D = N/V --> V = N/D --> L = (N/D)**(1/3)
boxSize = nPartDim/(density)**(1./3.)
end subroutine readData

end module read_data_module
