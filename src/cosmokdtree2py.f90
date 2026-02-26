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
    ! BUILDING FUNCTION (PARALLEL)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pybuild_kdtree(points, leaf, boxsize, c_ptr_tree)
        use omp_lib
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
    ! KNN SEARCH IN PARALLEL
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pyknn_search(c_ptr_tree, targett, ndim, ntar, k, sorted, dist, idx)
        use omp_lib
        implicit none
        !in
        integer*8, intent(in) :: c_ptr_tree
        integer, intent(in) :: ntar, k, ndim
        real*8, intent(in) :: targett(ntar,ndim)
        logical, intent(in) :: sorted
        !out
        real*8, intent(out) :: dist(ntar,k)
        integer*8, intent(out) :: idx(ntar,k)
        !local
        integer :: i
        type(KDTreeResult) :: query
        type(KDTreeNode), pointer :: tree
        type(c_ptr) :: temp_c_ptr

        !f2py intent(in) :: ndim, ntar, k
        !f2py intent(in) :: c_ptr_tree, targett, sorted
        !f2py intent(out) :: dist, idx
        !f2py depend(ntar,k) :: dist, idx
        !f2py depend(ntar,ndim) :: targett

        !to c pointer from integer representation
        temp_c_ptr = transfer(c_ptr_tree, temp_c_ptr)

        !to fortran pointer from c pointer
        call c_f_pointer(temp_c_ptr, tree)

        !$omp parallel shared(tree, targett, k, dist, idx) private(i, query)
        !$omp do
        do i = 1, ntar
             query = knn_search(tree, targett(i,:), k, sorted)
             dist(i,:) = query%dist
             idx(i,:) = query%idx
        enddo
        !$omp end do
        !$omp end parallel

    end subroutine pyknn_search
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !WITH BALL_SEARCH IS WAY MORE COMPLEX, ALLOCATABLE RESULTS NEED A 2-STAGE CALL 
    !QUERY RESULT IS TREATED AS AN OPAQUE POINTER SUCH AS TREE TO PASS TO PYTHON
    !TO ENSURE EFFICIENCY, ONLY 1 CALL TO THE FUNCTION IS MADE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! BALL SEARCH
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pyball_search_call_1(c_ptr_tree, targett, ndim, ntar, radius, sorted, &
                                    nball, nballmax, c_ptr_query)
        use omp_lib 
        implicit none
        !in
        integer, intent(in) :: ndim, ntar
        integer*8, intent(in) :: c_ptr_tree
        real*8, intent(in) :: targett(ntar,ndim)
        real*8, intent(in) :: radius
        logical, intent(in) :: sorted
        !out
        integer :: i
        integer, intent(out) :: nballmax, nball(ntar)
        integer*8, intent(out) :: c_ptr_query(ntar)
        !local
        type(KDTreeResult), pointer :: query
        type(KDTreeNode), pointer :: tree
        type(c_ptr) :: temp_c_ptr

        !f2py intent(in) :: ndim, ntar
        !f2py intent(in) :: c_ptr_tree, targett, radius, sorted
        !f2py intent(out) :: nball, nballmax, c_ptr_query
        !f2py depend(ntar) :: nball, c_ptr_query
        !f2py depend(ntar,ndim) :: targett

        !to c pointer from integer representation
        temp_c_ptr = transfer(c_ptr_tree, temp_c_ptr)

        !to fortran pointer from c pointer
        call c_f_pointer(temp_c_ptr, tree)

        !$omp parallel shared(ntar, tree, targett, radius, sorted, nball, c_ptr_query) &
        !$omp private(i, query, temp_c_ptr)
        !$omp do
        do i = 1, ntar
            allocate(query)
            query = ball_search(tree, targett(i,:), radius, sorted)
            nball(i) = size(query%idx)

            !opaque pointer to pass to python
            temp_c_ptr = c_loc(query)
            c_ptr_query(i) = transfer(temp_c_ptr, c_ptr_query(i))
        enddo
        !$omp end do
        !$omp end parallel

        nballmax = maxval(nball)

    end subroutine pyball_search_call_1

    subroutine pyball_search_call_2(ndim, ntar, c_ptr_query, nballmax, nball, dist, idx)
        implicit none
        !in
        integer*8, intent(in) :: c_ptr_query(ntar)
        integer, intent(in) :: nball(ntar)
        integer, intent(in) :: ndim, ntar, nballmax
        !out
        real*8, intent(out) :: dist(ntar, nballmax)
        integer*8, intent(out) :: idx(ntar, nballmax)
        !local
        integer :: i
        type(KDTreeResult), pointer :: query
        type(c_ptr) :: temp_c_ptr

        !f2py intent(in) :: c_ptr_query, nball, ndim, ntar, nballmax
        !f2py intent(out) :: dist, idx
        !f2py depend(ntar, nballmax) :: dist, idx
        !f2py depend(ntar) :: nball, c_ptr_query

        !$omp parallel shared(c_ptr_query, nball, dist, idx) private(i, query, temp_c_ptr)
        !$omp do
        do i = 1, ntar
            if (nball(i) > 0) then
                !to c pointer from integer representation
                temp_c_ptr = transfer(c_ptr_query(i), temp_c_ptr)

                !to fortran pointer from c pointer
                call c_f_pointer(temp_c_ptr, query)

                dist(i,1:nball(i)) = query%dist
                idx(i,1:nball(i)) = query%idx

                !deallocate query
                deallocate(query)
            else
                dist(i,:) = 0.
                idx(i,:) = 0
            endif
        enddo
        !$omp end do
        !$omp end parallel
        
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
