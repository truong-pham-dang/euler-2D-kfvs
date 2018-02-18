subroutine solver
    implicit none
    
    call solver_kfvs 
end subroutine solver
    
subroutine solver_kfvs
    use mod_solver
    use mod_struct_to_array
    use mod_write_vtk, only: write_solution_vtk
    implicit none 
    integer :: n
    
    call struct_to_array 
    
    call donnee_initiale
    call allocate_vardummy
    call conditions_aux_limites
    call assign_lr_cell 
    call calcul_derived_quantities 
    call calcul_conservative_vector 
    
    call timestep 
    
    !--- temps de simulation et parametres de sorties
    tmax=0.3828823925d-2
    nmax=floor(tmax/dt)
    
    ! debug
    nmax = 500000
    
    !--- evolution
    do n=1,nmax
        
        !-- calcul des flux
        call calcul_flux 
        
        !-- iteration en temps
        call calcul_rhs
        call euler_time_iteration
        
        !-- mise a jour
        vect_u=vect_unew
        
        !-- calcul de rho,ux,uy,t
        call calcul_rho_ux_uy_t
        
        !-- mise a jour cl
        call conditions_aux_limites
        
        !-- mise à jour des quantités dérivées
        call calcul_derived_quantities
        
        !-- sauvegarde resultats (format vtk, lisible par Paraview)
        if (mod(n,10000) == 0) then
            write(*,*) 'Writing solution file at iteration ', n, '...'
            call write_solution_vtk(n)
        endif 
    enddo 
    
end subroutine solver_kfvs