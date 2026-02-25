!This is a wrapper using iso_c_binding to
!easily pass "opaque" derived type data from fortran to python
!as this data is only really accessed inside 

!Double-layer used to handle derived type data and optional inputs
!--------------------------------------------------------------
!cosmokdtree <-> cosmokdtree2py <-> python_wrapper <-> python
!--------------------------------------------------------------

!REAL*8 precision and INT*8 integers are used for better compatibility with python's
!  float64 and int64 types

module pycosmokdtree
    use cosmokdtree
    use iso_c_binding
    implicit none

    public :: pybuild_kdtree, pyknn_search, pyball_search_call_1, pyball_search_call_2, &
                pybox_search_call_1, pybox_search_call_2

contains
    !THIS MODULE ONLY CONTAINS WRAPPER FUNCTIONS FOR F2PY
    !actual kdtree is in cosmokdtree.f90 file


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! BUILDING FUNCTION
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pybuild_kdtree(points, leaf, boxsize, c_ptr_tree)
        implicit none
        !in
        real*8, intent(in) :: points(:,:)
        real*8, intent(in) :: boxsize(:)
        integer, intent(in) :: leaf
        !local
        type(KDTreeNode), pointer :: tree
        type(c_ptr) :: temp_c_ptr
        logical :: periodic_in
        integer :: i
        !out
        integer*8, intent(out) :: c_ptr_tree

        !f2py intent(out) :: c_ptr_tree
        !f2py intent(in) :: points, leaf, boxsize

        !check boxsize and periodicity
#if periodic == 1
        periodic_in = .false.
        do i = 1, size(boxsize)
            if ( boxsize(i) /= -1. ) then
                periodic_in = .true.
                exit
            endif
        enddo

        if ( periodic_in .eqv. .false. ) then
            STOP "Error: provide boxsize, periodic boundary conditions are enabled."
        endif
#endif

        !building function
#if periodic == 1
        tree => build_kdtree(points, boxsize, leaf)
#else
        tree => build_kdtree(points, leaf)
#endif
        
        !C adress pointing to the fortran complexx data
        temp_c_ptr = c_loc(tree)

        !to integer representation of adress to f2py it
        c_ptr_tree = transfer(temp_c_ptr, c_ptr_tree)
    end subroutine pybuild_kdtree
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DEALLOCATE TREE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pydeallocate_tree(c_ptr_tree)
        implicit none
        !in
        integer*8, intent(in) :: c_ptr_tree
        !local
        type(KDTreeNode), pointer :: tree
        type(c_ptr) :: temp_c_ptr

        !f2py intent(in) :: c_ptr_tree

        !to c pointer from integer representation
        temp_c_ptr = transfer(c_ptr_tree, temp_c_ptr)

        !to fortran pointer from c pointer
        call c_f_pointer(temp_c_ptr, tree)

        call deallocate_kdtree(tree)
    end subroutine pydeallocate_tree
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! KNN SEARCH
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pyknn_search(c_ptr_tree, targett, k, sorted, dist, idx)
        implicit none
        !in
        integer*8, intent(in) :: c_ptr_tree
        real*8, intent(in) :: targett(:)
        integer, intent(in) :: k
        logical, intent(in) :: sorted
        !out
        real*8, intent(out) :: dist(k)
        integer*8, intent(out) :: idx(k)
        !local
        type(KDTreeResult) :: query
        type(KDTreeNode), pointer :: tree
        type(c_ptr) :: temp_c_ptr

        !f2py intent(in) :: c_ptr_tree, targett, k, sorted
        !f2py intent(out) :: dist, idx
        !f2py depend(k) :: dist, idx

        !to c pointer from integer representation
        temp_c_ptr = transfer(c_ptr_tree, temp_c_ptr)

        !to fortran pointer from c pointer
        call c_f_pointer(temp_c_ptr, tree)

        query = knn_search(tree, targett, k, sorted)

        dist = query%dist
        idx = query%idx
    end subroutine pyknn_search
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !WITH BALL_SEARCH IS WAY MORE COMPLEX, ALLOCATABLE RESULTS NEED A 2-STAGE CALL 
    !QUERY RESULT IS TREATED AS AN OPAQUE POINTER SUCH AS TREE TO PASS TO PYTHON
    !TO ENSURE EFFICIENCY, ONLY 1 CALL TO THE FUNCTION IS MADE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! BALL SEARCH
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pyball_search_call_1(c_ptr_tree, targett, radius, sorted, nball, c_ptr_query)
        implicit none
        !in
        integer*8, intent(in) :: c_ptr_tree
        real*8, intent(in) :: targett(:)
        real*8, intent(in) :: radius
        logical, intent(in) :: sorted
        !out
        integer, intent(out) :: nball
        integer*8, intent(out) :: c_ptr_query
        !local
        type(KDTreeResult), pointer :: query
        type(KDTreeNode), pointer :: tree
        type(c_ptr) :: temp_c_ptr

        !f2py intent(in) :: c_ptr_tree, target, radius, sorted
        !f2py intent(out) :: nball, c_ptr_query

        !to c pointer from integer representation
        temp_c_ptr = transfer(c_ptr_tree, temp_c_ptr)

        !to fortran pointer from c pointer
        call c_f_pointer(temp_c_ptr, tree)

        allocate(query)
        query = ball_search(tree, targett, radius, sorted)
    
        nball = size(query%idx)

        !opaque pointer to pass to python
        temp_c_ptr = c_loc(query)
        c_ptr_query = transfer(temp_c_ptr, c_ptr_query)

    end subroutine pyball_search_call_1

    subroutine pyball_search_call_2(c_ptr_query, nball, dist, idx)
        implicit none
        !in
        integer*8, intent(in) :: c_ptr_query
        integer, intent(in) :: nball
        !out
        real*8, intent(out) :: dist(nball)
        integer*8, intent(out) :: idx(nball)
        !local
        type(KDTreeResult), pointer :: query
        type(c_ptr) :: temp_c_ptr

        !f2py intent(in) :: c_ptr_query, nball
        !f2py intent(out) :: dist, idx
        !f2py depend(nball) :: dist, idx

        !to c pointer from integer representation
        temp_c_ptr = transfer(c_ptr_query, temp_c_ptr)
        
        !to fortran pointer from c pointer
        call c_f_pointer(temp_c_ptr, query)

        dist = query%dist
        idx = query%idx

        !deallocate query
        deallocate(query)
        
    end subroutine pyball_search_call_2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! BOX SEARCH
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pybox_search_call_1(c_ptr_tree, box, nbox, c_ptr_query)
        implicit none
        !in
        integer*8, intent(in) :: c_ptr_tree
        real*8, intent(in) :: box(:)
        !out
        integer, intent(out) :: nbox
        integer*8, intent(out) :: c_ptr_query
        !local
        type(KDTreeResult), pointer :: query
        type(KDTreeNode), pointer :: tree
        type(c_ptr) :: temp_c_ptr

        !f2py intent(in) :: c_ptr_tree, box
        !f2py intent(out) :: nbox, c_ptr_query

        !to c pointer from integer representation
        temp_c_ptr = transfer(c_ptr_tree, temp_c_ptr)   

        !to fortran pointer from c pointer
        call c_f_pointer(temp_c_ptr, tree)  

        !check dimensionality against box
        if (size(box) /= 2*size(tree%maxbounds)) then
            STOP "Error: box size does not match tree dimensionality."
        endif

        allocate(query)
        query = box_search(tree, box)   

        nbox = size(query%idx)

        !opaque pointer to pass to python
        temp_c_ptr = c_loc(query)
        c_ptr_query = transfer(temp_c_ptr, c_ptr_query)
    end subroutine pybox_search_call_1

    subroutine pybox_search_call_2(c_ptr_query, nbox, idx)
        implicit none
        !in
        integer*8, intent(in) :: c_ptr_query
        integer, intent(in) :: nbox
        !out
        integer*8, intent(out) :: idx(nbox)
        !local
        type(KDTreeResult), pointer :: query
        type(c_ptr) :: temp_c_ptr

        !f2py intent(in) :: c_ptr_query, nbox
        !f2py intent(out) :: idx
        !f2py depend(nbox) :: idx

        !to c pointer from integer representation
        temp_c_ptr = transfer(c_ptr_query, temp_c_ptr)

        !to fortran pointer from c pointer
        call c_f_pointer(temp_c_ptr, query)

        idx = query%idx

        !deallocate query
        deallocate(query)

    end subroutine pybox_search_call_2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module pycosmokdtree
