
subroutine pre_processing
    use mod_read_gmsh
    use mod_write_vtk
    use mod_cell_2D
    use mod_detect_nearest_neighbor
    use mod_fvm_face_2D
    implicit none


    call read_mesh
    call construct_id_nodes 
    call construct_cells
    call detect_neighbor
    call write_mesh_vtk
    call write_mesh_tecplot    
    call calcul_vol_cells    
    call calcul_area_cent_faces    
    call assign_id_face_v2   
    call allocate_face_2D  
    call assign_face_2D 

end subroutine pre_processing 

