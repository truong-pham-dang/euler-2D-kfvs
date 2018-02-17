module mod_struct_to_array
    use mod_cell_2D
    use mod_fvm_face_2D
    implicit none
    
    real(8), dimension(:), allocatable   :: vol, norm_x, norm_y, len_norm
    integer, dimension(:), allocatable   :: bc_typ
    integer, dimension(:,:), allocatable :: lr_cell
    
    contains
!----------------------------------------------------------------------
    subroutine struct_to_array
    implicit none
    
    integer :: icel, ifac
    
    if (.not. allocated(vol)) then
        allocate(vol(nbelm))
    endif 
    
    if (.not. allocated(norm_x)) then
        allocate(norm_x(nbfaces))
    endif 
    
    if (.not. allocated(norm_y)) then
        allocate(norm_y(nbfaces))
    endif
    
    if (.not. allocated(len_norm)) then
        allocate(len_norm(nbfaces))
    endif
    
    if (.not. allocated(bc_typ)) then
        allocate(bc_typ(nbfaces))
    endif
    
    if (.not. allocated(lr_cell)) then
        allocate(lr_cell(nbfaces,2))
    endif
    
    do icel = 1, nbelm
        vol(icel) = cell(icel)%p%vol
    enddo 
    
    do ifac = 1, nbfaces
        norm_x(ifac)   = face_2D(ifac)%f%normal_vector%x
        norm_y(ifac)   = face_2D(ifac)%f%normal_vector%y
        len_norm(ifac) = face_2D(ifac)%f%len_nor
        bc_typ(ifac)   = face_2D(ifac)%f%bc_typ 
    enddo 
    
    end subroutine struct_to_array
end module