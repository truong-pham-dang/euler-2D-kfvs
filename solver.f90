subroutine solver
    implicit none
    
    call solver_kfvs 
end subroutine solver
    
subroutine solver_kfvs
    use mod_solver
    implicit none 
    
    call donnee_initiale
    call allocate_vardummy
    call conditions_aux_limites
    call calcul_derived_quantities 
    call calcul_conservative_vector 
    
    call timestep 
    
end subroutine solver_kfvs