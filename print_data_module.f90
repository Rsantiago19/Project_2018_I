module print_data_module
implicit none

contains

subroutine print_data(time, V, KE, Tins, mom, un)
implicit none
integer, intent(in)				:: un
real(8), intent(in)				:: time, V, KE, Tins, mom
! Imprimeix en un fitxer que te per unitat 'un' el temps, la energia potencial,
! la energia cinètica, la energia mecànica, el moment total i la temperatura
write(un,*) time, V, KE, KE + V, mom, Tins
end subroutine print_data

end module print_data_module
