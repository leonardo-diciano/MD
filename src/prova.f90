module prova
use force_field_mod

contains

subroutine provaaa(positions,tot_pot,forces)
    real, intent(inout) :: tot_pot
    real, allocatable, intent(inout) :: forces(:,:), positions(:,:)

    positions=positions*1.05
    
    CALL get_energy_gradient(positions,tot_pot,forces)

end subroutine

end module