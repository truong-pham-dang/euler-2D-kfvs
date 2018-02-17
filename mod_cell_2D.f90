module mod_cell_2D
    use mod_read_gmsh
    use mod_point_util
    implicit none

    type face
        type(point)            :: p1,p2
        integer(4)             :: bc_typ = 0
        integer(4)             :: idface = 0
        real(8)                :: area = 0.0d0
        type(point)            :: centroid
    end type face

    type cell_2D
        integer(4)  :: ident = 0
        type(point) :: vertex(4)

        type(face),dimension(4) :: faces

        type(cell_2D),pointer :: neighbor1 => null()
        type(cell_2D),pointer :: neighbor2 => null()
        type(cell_2D),pointer :: neighbor3 => null()
        type(cell_2D),pointer :: neighbor4 => null()
        
        real(8)     :: vol = 0.0d0

    end type cell_2D
     
    type list_cell_2D
        type(cell_2D),pointer :: p
    end type

    type(list_cell_2D),dimension(:),allocatable :: cell
    
    integer(4)      :: nbfaces
    
    private :: sort_vertex

    contains
!----------------------------------------------------------------------    
       subroutine construct_cells
       implicit none

       integer(4) :: i,j
       integer(4) :: idnode
       
       integer(4) :: idnode1, idnode2
       
       type(cell_2D), pointer        :: pcell
       type(node_ident_msh), pointer :: pmsh

       allocate(cell(1:nbelm))

       do i = 1, nbelm
          allocate(cell(i)%p)
       end do

       ! Assign vertexes of a cell

       do i = 1, nbelm
         cell(i)%p%ident = i
         cell(i)%p%vertex(1)%ident = id_nodes(i)%pn%id_node(6)
         idnode = cell(i)%p%vertex(1)%ident
         do j = 1, nbnode
            if (idnode == coord_nodes(j)%p%ident) then
                cell(i)%p%vertex(1)%x = coord_nodes(j)%p%x
                cell(i)%p%vertex(1)%y = coord_nodes(j)%p%y
            end if
         end do

         cell(i)%p%vertex(2)%ident = id_nodes(i)%pn%id_node(7)
         idnode = cell(i)%p%vertex(2)%ident
         do j = 1, nbnode
            if (idnode == coord_nodes(j)%p%ident) then
                cell(i)%p%vertex(2)%x = coord_nodes(j)%p%x
                cell(i)%p%vertex(2)%y = coord_nodes(j)%p%y
            end if
         end do

         cell(i)%p%vertex(3)%ident = id_nodes(i)%pn%id_node(8)
         idnode = cell(i)%p%vertex(3)%ident
         do j = 1, nbnode
            if (idnode == coord_nodes(j)%p%ident) then
                cell(i)%p%vertex(3)%x = coord_nodes(j)%p%x
                cell(i)%p%vertex(3)%y = coord_nodes(j)%p%y
            end if
         end do

         cell(i)%p%vertex(4)%ident = id_nodes(i)%pn%id_node(9)
         idnode = cell(i)%p%vertex(4)%ident
         do j = 1, nbnode
            if (idnode == coord_nodes(j)%p%ident) then
                cell(i)%p%vertex(4)%x = coord_nodes(j)%p%x
                cell(i)%p%vertex(4)%y = coord_nodes(j)%p%y
            end if
         end do

       end do
       
       ! Sort vertexes for each cell
       ! sort_vertex failed
       !do i = 1, nbelm
       !    call sort_vertex(cell(i)%p%vertex)
       !enddo 

       ! Assign faces for each cell

       do i = 1, nbelm
         cell(i)%p%faces(1)%p1%x = cell(i)%p%vertex(1)%x
         cell(i)%p%faces(1)%p1%y = cell(i)%p%vertex(1)%y

         cell(i)%p%faces(1)%p2%x = cell(i)%p%vertex(2)%x
         cell(i)%p%faces(1)%p2%y = cell(i)%p%vertex(2)%y

         cell(i)%p%faces(2)%p1%x = cell(i)%p%vertex(2)%x
         cell(i)%p%faces(2)%p1%y = cell(i)%p%vertex(2)%y

         cell(i)%p%faces(2)%p2%x = cell(i)%p%vertex(3)%x
         cell(i)%p%faces(2)%p2%y = cell(i)%p%vertex(3)%y

         cell(i)%p%faces(3)%p1%x = cell(i)%p%vertex(3)%x
         cell(i)%p%faces(3)%p1%y = cell(i)%p%vertex(3)%y

         cell(i)%p%faces(3)%p2%x = cell(i)%p%vertex(4)%x
         cell(i)%p%faces(3)%p2%y = cell(i)%p%vertex(4)%y

         cell(i)%p%faces(4)%p1%x = cell(i)%p%vertex(4)%x
         cell(i)%p%faces(4)%p1%y = cell(i)%p%vertex(4)%y

         cell(i)%p%faces(4)%p2%x = cell(i)%p%vertex(1)%x
         cell(i)%p%faces(4)%p2%y = cell(i)%p%vertex(1)%y

       end do
       
       ! Assign boundary condition for each face of a cell 
       
       do i = 1, nbelm
           pcell => cell(i)%p 
           idnode1 = pcell%vertex(1)%ident 
           idnode2 = pcell%vertex(2)%ident 
           
           loop_msh1: do j = 1, nbel_msh
               pmsh => id_nodes_msh(j)%pn 
               if (pmsh%elem_typ == 27) then ! I write a solver for NACA 0012 HIOCFD firstly
                   if (idnode1 == pmsh%id_node(6) .and. idnode2 == pmsh%id_node(7)) then
                       pcell%faces(1)%bc_typ = pmsh%tag1
                       cycle loop_msh1
                   else if (idnode1 == pmsh%id_node(7) .and. idnode2 == pmsh%id_node(6)) then
                       pcell%faces(1)%bc_typ = pmsh%tag1
                       cycle loop_msh1
                   endif 
               endif
           enddo loop_msh1
           
           idnode1 = pcell%vertex(2)%ident 
           idnode2 = pcell%vertex(3)%ident 
           
           loop_msh2: do j = 1, nbel_msh
               pmsh => id_nodes_msh(j)%pn 
               if (pmsh%elem_typ == 27) then ! I write a solver for NACA 0012 HIOCFD firstly
                   if (idnode1 == pmsh%id_node(6) .and. idnode2 == pmsh%id_node(7)) then
                       pcell%faces(2)%bc_typ = pmsh%tag1
                       cycle loop_msh2
                   else if (idnode1 == pmsh%id_node(7) .and. idnode2 == pmsh%id_node(6)) then
                       pcell%faces(2)%bc_typ = pmsh%tag1
                       cycle loop_msh2
                   endif 
               endif
           enddo loop_msh2
           
           idnode1 = pcell%vertex(3)%ident 
           idnode2 = pcell%vertex(4)%ident 
           
           loop_msh3: do j = 1, nbel_msh
               pmsh => id_nodes_msh(j)%pn 
               if (pmsh%elem_typ == 27) then ! I write a solver for NACA 0012 HIOCFD firstly
                   if (idnode1 == pmsh%id_node(6) .and. idnode2 == pmsh%id_node(7)) then
                       pcell%faces(3)%bc_typ = pmsh%tag1
                       cycle loop_msh3
                   else if (idnode1 == pmsh%id_node(7) .and. idnode2 == pmsh%id_node(6)) then
                       pcell%faces(3)%bc_typ = pmsh%tag1
                       cycle loop_msh3
                   endif 
               endif
           enddo loop_msh3
           
           idnode1 = pcell%vertex(4)%ident 
           idnode2 = pcell%vertex(1)%ident 
           
           loop_msh4: do j = 1, nbel_msh
               pmsh => id_nodes_msh(j)%pn 
               if (pmsh%elem_typ == 27) then ! I write a solver for NACA 0012 HIOCFD firstly
                   if (idnode1 == pmsh%id_node(6) .and. idnode2 == pmsh%id_node(7)) then
                       pcell%faces(4)%bc_typ = pmsh%tag1
                       cycle loop_msh4
                   else if (idnode1 == pmsh%id_node(7) .and. idnode2 == pmsh%id_node(6)) then
                       pcell%faces(4)%bc_typ = pmsh%tag1
                       cycle loop_msh4
                   endif 
               endif
           enddo loop_msh4
           
       enddo
       
       end subroutine
!----------------------------------------------------------------------  
       subroutine calcul_vol_cells
       implicit none 
       integer                :: i
       type(cell_2D), pointer :: pc
       real(8)                :: vol, x1, x2, x3, x4, y1, y2, y3, y4
       
       do i = 1, nbelm
           pc => cell(i)%p
           
           x1     = pc%vertex(1)%x
           x2     = pc%vertex(2)%x
           x3     = pc%vertex(3)%x
           x4     = pc%vertex(4)%x
           
           y1     = pc%vertex(1)%y
           y2     = pc%vertex(2)%y
           y3     = pc%vertex(3)%y
           y4     = pc%vertex(4)%y
           
           vol    = 0.5d0 * ( (x1 - x3)*(y2 - y4) + (x4 - x2)*(y1 - y3) )
           pc%vol = vol
       enddo 
       end subroutine calcul_vol_cells
!----------------------------------------------------------------------
       subroutine calcul_area_cent_faces
       use mod_vector_algebra
       implicit none 
       integer                :: i, j
       type(cell_2D), pointer :: pcel
       type(face),    pointer :: pfac
       type(vector_t)         :: vec
       
       do i = 1, nbelm
           pcel => cell(i)%p
           do j = 1, 4
               pfac => pcel%faces(j)
               vec           = vector_t( pfac%p1%x - pfac%p2%x, pfac%p1%y - pfac%p2%y, 0.0d0 )
               pfac%area     = .abs.(vec)
               pfac%centroid%x =  0.5d0 * (pfac%p1%x + pfac%p2%x)
               pfac%centroid%y =  0.5d0 * (pfac%p1%y + pfac%p2%y)
           enddo 
       enddo 
       end subroutine calcul_area_cent_faces
!----------------------------------------------------------------------  
       subroutine assign_id_face
       use mod_point_util
       implicit none
       integer                :: i, j, k
       type(cell_2D), pointer :: pc, pc2
       type(face),    pointer :: pfac, pfac2
       
       nbfaces = 0
       do i = 1, nbelm
           pc => cell(i)%p
               pfac => pc%faces(1)
               if (pfac%idface == 0) then ! a face that has not been indexed yet
               if (associated(pc%neighbor1)) then
                  pc2 => pc%neighbor1
                  if (associated(pc%neighbor1%neighbor3, pc)) then 
                    if (pfac%bc_typ == 0) then ! normal face
                        pfac2 => pc2%faces(3)
                        if ((pfac%centroid .eq. pfac2%centroid) .and. abs(pfac%area - pfac2%area) <= 1.0d-9) then
                            if (pfac2%idface == 0) then
                                nbfaces = nbfaces + 1
                                pfac%idface = nbfaces
                                pfac2%idface = nbfaces
                            endif
                        endif
                    endif
                  endif 
                else if (.not. associated(pc%neighbor1)) then 
                   if (pfac%bc_typ /= 0) then ! boundary face
                        nbfaces = nbfaces + 1
                        pfac%idface = nbfaces 
                   endif
                endif
               endif 
               
               pfac => pc%faces(2)
               if (pfac%idface == 0) then ! a face that has not been indexed yet
               if (associated(pc%neighbor2)) then
                  pc2 => pc%neighbor2
                  if (associated(pc%neighbor2%neighbor4, pc)) then 
                    if (pfac%bc_typ == 0) then ! normal face
                        pfac2 => pc2%faces(4)
                        if ((pfac%centroid .eq. pfac2%centroid) .and. abs(pfac%area - pfac2%area) <= 1.0d-9) then
                            if (pfac2%idface == 0) then
                                nbfaces = nbfaces + 1
                                pfac%idface = nbfaces
                                pfac2%idface = nbfaces
                            endif
                        endif
                    endif
                  endif 
                else if (.not. associated(pc%neighbor2)) then 
                   if (pfac%bc_typ /= 0) then ! boundary face
                        nbfaces = nbfaces + 1
                        pfac%idface = nbfaces 
                   endif
                endif
               endif
               
               pfac => pc%faces(3)
               if (pfac%idface == 0) then ! a face that has not been indexed yet
               if (associated(pc%neighbor3)) then
                  pc2 => pc%neighbor3
                  if (associated(pc%neighbor3%neighbor1, pc)) then 
                    if (pfac%bc_typ == 0) then ! normal face
                        pfac2 => pc2%faces(1)
                        if ((pfac%centroid .eq. pfac2%centroid) .and. abs(pfac%area - pfac2%area) <= 1.0d-9) then
                            if (pfac2%idface == 0) then
                                nbfaces = nbfaces + 1
                                pfac%idface = nbfaces
                                pfac2%idface = nbfaces
                            endif
                        endif
                    endif
                  endif 
                else if (.not. associated(pc%neighbor3)) then 
                   if (pfac%bc_typ /= 0) then ! boundary face
                        nbfaces = nbfaces + 1
                        pfac%idface = nbfaces 
                   endif
                endif
               endif
               
               pfac => pc%faces(4)
               if (pfac%idface == 0) then ! a face that has not been indexed yet
               if (associated(pc%neighbor4)) then
                  pc2 => pc%neighbor4
                  if (associated(pc%neighbor4%neighbor2, pc)) then 
                    if (pfac%bc_typ == 0) then ! normal face
                        pfac2 => pc2%faces(2)
                        if ((pfac%centroid .eq. pfac2%centroid) .and. abs(pfac%area - pfac2%area) <= 1.0d-9) then
                            if (pfac2%idface == 0) then
                                nbfaces = nbfaces + 1
                                pfac%idface = nbfaces
                                pfac2%idface = nbfaces
                            endif
                        endif
                    endif
                  endif 
                else if (.not. associated(pc%neighbor4)) then 
                   if (pfac%bc_typ /= 0) then ! boundary face
                        nbfaces = nbfaces + 1
                        pfac%idface = nbfaces 
                   endif
                endif
                endif
       enddo 
       end subroutine assign_id_face
!----------------------------------------------------------------------  
       subroutine sort_vertex(vertices)
       !Reference: http://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/sorting.f90                  
       implicit none 
       
       type(point) :: vertices(:)
       
       integer     :: i, location
       
       do i = 1, size(vertices) - 1
           location = findminimum(vertices(i),i,size(vertices))
           call swap_point(vertices(i),vertices(location))
       enddo 
       
       if (vertices(1)%x > vertices(2)%x) then
           call swap_point(vertices(1),vertices(2))
       endif 
       
       if (vertices(4)%x > vertices(3)%x) then
           call swap_point(vertices(4),vertices(3))
       endif
       
       contains
          integer function findminimum(in_point,start,endp)
          implicit none
          type(point)         :: in_point
          integer, intent(in) :: start, endp
          
          integer             :: location, i
          type(point)         :: minimum
          
          minimum = vertices(start)
          location = start
          do i = start + 1, endp
              if (vertices(i)%y < minimum%y) then
                  minimum = vertices(i)
                  location = i
              endif 
          enddo 
          findminimum = location
          
          end function findminimum
       end subroutine sort_vertex
!----------------------------------------------------------------------
end module
