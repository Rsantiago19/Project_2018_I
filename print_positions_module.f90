module print_positions_module
implicit none
contains

subroutine print_positions(un, nPart, pos, time)
implicit none
integer, intent(in)                             :: un, nPart
real(8), dimension(nPart,3), intent(in)         :: pos
real(8), intent(in)                             :: time
integer                                         :: i

write(un,*) nPart
write(un,*) 'Simulation Time = ', time,'seconds'
do i = 1, nPart, 1
        write(un,*) 'C', pos(i,:)
end do

end subroutine print_positions

end module print_positions_module
