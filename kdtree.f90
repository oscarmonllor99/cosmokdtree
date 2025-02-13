!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! kd-tree module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! GRUP DE COSMOLOGIA COMPUTACIONAL (GCC) UNIVERSITAT DE VALÈNCIA
! Author: Óscar Monllor Berbegal
! Date: 30/01/2025
! Last update: 13/02/2025
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Brief description:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! - This module implements an optimized parallel kd-tree construction and 
!   search for k-nearest neighbours and points within a given radius.
!
! - Quickselect is used to find the point that splits the space
!   in two halves along a given axis (median). Use of median ensures 
!   that the tree is balanced.
!
! - The tree is built recursively, with the splitting axis changing
!   at each level (x, y, z, x, y, z, ...). The tree is built in parallel 
!   using OpenMP tasks.
!
! - The search for k-nearest neighbours uses an insertion shiftdown 
!   too quickly sort and replace the nearest neighbours found.
!
! - quicksort is used to sort the distances and indices of    
!   points within a given radius
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Pending improvements:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! - Quicksort may not be the most efficient way to sort the distances
!   and indices of points within a given radius. They could be sorted
!   on the fly during the search by shifting the elements, just like
!   in the k-nearest neighbour search.
!
! - 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!#######################################################
module kdtree_mod
!#######################################################
    implicit none
    private  
    public :: build_kdtree_init, KDTreeNode, KDTreeResult, knn_search_init, ball_search_init
    !+++++++++++++++++++++++++++++++
    !++++ Type definitions
    !+++++++++++++++++++++++++++++++
    type :: KDTreeNode
        !basic ---------------------
        real :: point(3)       ! 3D point (x, y, z)
        integer(kind=8) :: index       ! Index of the point within the original array
        integer :: axis        ! Splitting axis (0 for x, 1 for y, 2 for z)
        type(KDTreeNode), pointer :: left => null()  ! Left child
        type(KDTreeNode), pointer :: right => null() ! Right child
        !leaf, for faster search and building
        integer :: is_leaf     ! Flag to indicate if the node is a leaf
        real, pointer :: leaf_points(:, :) => null()  ! Points in the leaf (for leaf nodes)
        integer(kind=8), pointer :: leaf_indices(:) => null()  ! Indices of points in the leaf
    end type KDTreeNode

    type :: KDTreeResult
        integer(kind=8), allocatable :: idx(:)
        real, allocatable :: dist(:)
    end type KDTreeResult
    !+++++++++++++++++++++++++++++++

contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Initialize kd-tree construction
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function build_kdtree_init(x, y, z) result(tree)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      use omp_lib
      implicit none
      real, intent(in) :: x(:), y(:), z(:)
      real, allocatable :: points(:, :)
      integer(kind=8), allocatable :: indices(:)
      integer(kind=8) :: n, i
      integer :: depth, max_depth, nproc, leafsize
      type(KDTreeNode), pointer :: tree
      
      ! Enable nested parallelism
      call omp_set_nested(.true.) 
    
      ! Number of points
      n = size(x, kind=8)

      ! 3D points array
      allocate(points(n, 3))
      allocate(indices(n))

      points(:, 1) = x
      points(:, 2) = y
      points(:, 3) = z

      ! Initialize global indices
      indices = [(i, i=1, n)]
        
      ! Init depth = 0
      depth = 0

      !$OMP PARALLEL
      !$OMP SINGLE
      nproc = omp_get_num_threads()
      !$OMP END SINGLE
      !$OMP END PARALLEL

      ! Build KD-tree
      max_depth = 2 + compute_max_depth(omp_get_max_threads())
      write(*,*) "Parallel max depth:", max_depth

      ! Leafsize scaling with the number of points
      leafsize = int(real(n)**0.333 / 4.)
      leafsize = max(leafsize, 1)

      tree => build_kdtree(points, indices, depth, max_depth, leafsize)

      deallocate(points, indices)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end function build_kdtree_init
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !Compute maximum depth of the tree for parallelism (Ncores)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function compute_max_depth(nproc) result(max_depth)
        implicit none
        integer, intent(in) :: nproc
        integer :: max_depth
        !nproc: number of processes/threads available.
        !max_depth Maximum depth which guarantees that there is at least one idle process.
        max_depth = int (log(dble(nproc)+0.1) / log(2.))
    end function compute_max_depth
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    recursive function build_kdtree(points,indices,depth,max_depth,leafsize) result(node)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    use omp_lib
    implicit none
    real, intent(inout) :: points(:, :)      ! 2D array of points
    integer(kind=8), intent(inout) :: indices(:) ! 1D array of indices
    integer, intent(in) :: depth         ! Current depth in the tree
    integer, intent(in) :: max_depth
    integer, intent(in) :: leafsize
    integer :: axis
    type(KDTreeNode), pointer :: node   ! New node to be created
    real :: kth_point(size(points, 2))
    integer(kind=8) :: median
    integer(kind=8) :: kth_index

    if (size(points, 1) == 0) then
        node => null()
        return
    end if

    !alternate axis across depth
    axis = mod(depth, 3)

    ! Find median and partition points
    median = size(points, 1, kind=8) / 2 + 1
    call quickselect(points, indices, median, axis, kth_point, kth_index)

    ! Allocate node
    allocate(node)
    node%point = kth_point
    node%index = kth_index
    node%axis = axis

    ! Check if this is a leaf node
    if (size(points, 1, kind=8) <= leafsize) then
        node%is_leaf = 1
        ! Store all points and indices in the leaf node
        allocate(node%leaf_points(size(points, 1, kind=8), size(points, 2)))
        allocate(node%leaf_indices(size(points, 1, kind=8)))
        node%leaf_points = points
        node%leaf_indices = indices
        node%left => null()
        node%right => null()
        return
    else
        node%is_leaf = 0
        node%point = kth_point
        node%index = kth_index
    end if

    ! Subtree construction (parallel at the top levels)
    if (depth < max_depth) then
    !$OMP PARALLEL num_threads(2**(max_depth-depth))
    !$OMP SINGLE

    !$OMP TASK
    node%left => build_kdtree(points(1:median-1,:),indices(1:median-1),depth+1,max_depth,leafsize)
    !$OMP END TASK

    !$OMP TASK
    node%right => build_kdtree(points(median+1:,:),indices(median+1:),depth+1,max_depth,leafsize)
    !$OMP END TASK

    ! Wait for the tasks to complete
    !$OMP TASKWAIT
    !$OMP END SINGLE
    !$OMP END PARALLEL
    else
    node%left => build_kdtree(points(1:median-1,:),indices(1:median-1),depth+1,max_depth,leafsize)
    node%right => build_kdtree(points(median+1:,:),indices(median+1:),depth+1,max_depth,leafsize)
    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end function build_kdtree
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Quickselect function to find the k-th smallest element along a specified axis
    ! Ensures all points below k are less than or equal to the k-th point
    ! and all points above k are greater than or equal to the k-th point with the
    ! specified axis.
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine quickselect(points, indices, k, axis, kth_point, kth_index)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        implicit none
        real, intent(inout) :: points(:, :)      ! 2D array of points (slice of the full array)
        integer(kind=8), intent(inout) :: indices(:) ! 1D array of indices (slice of the full array)
        integer(kind=8), intent(in) :: k         ! k-th smallest element to find
        integer, intent(in) :: axis             ! Axis to sort along (0 for x, 1 for y, 2 for z)
        real, intent(out) :: kth_point(size(points, 2))  ! The k-th smallest point
        integer(kind=8), intent(out) :: kth_index        ! Index of the k-th smallest point (within the original array)
        integer(kind=8) :: left, right, pivot_index

        left = 1
        right = size(points, 1, kind=8)
        do while (left <= right)
            ! Partition the array and get the pivot index
            pivot_index = partition(points, indices, left, right, axis)

            if (pivot_index == k) then
                ! Found the k-th smallest element
                kth_point = points(pivot_index, :)
                kth_index = indices(pivot_index)
                return
            else if (pivot_index < k) then
                ! Search the right subarray
                left = pivot_index + 1
            else
                ! Search the left subarray
                right = pivot_index - 1
            end if
        end do

        ! If the loop ends, return the k-th element
        kth_point = points(k, :)
        kth_index = indices(k)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine quickselect
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    ! Partition function for Quickselect
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function partition(points, indices, left, right, axis) result(pivot_index)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        implicit none
        real, intent(inout) :: points(:, :)      ! 2D array of points (slice of the full array)
        integer(kind=8), intent(inout) :: indices(:) ! 1D array of indices (slice of the full array)
        integer(kind=8), intent(in) :: left, right  ! Left and right bounds of the partition
        integer, intent(in) :: axis             ! Axis to sort along (0 for x, 1 for y, 2 for z)
        integer(kind=8) :: pivot_index, i, j
        real :: pivot_value, temp_point(size(points, 2))
        integer(kind=8) :: temp_index

        ! Use middle element as pivot
        pivot_value = points((left + right) / 2, axis + 1)

        i = left - 1

        do j = left, right - 1
            if (points(j, axis + 1) <= pivot_value) then
                i = i + 1
                ! Swap points(i, :) and points(j, :)
                temp_point = points(i, :)
                points(i, :) = points(j, :)
                points(j, :) = temp_point

                ! Swap indices(i) and indices(j)
                temp_index = indices(i)
                indices(i) = indices(j)
                indices(j) = temp_index
            end if
        end do

        ! Swap points(i+1, :) and points(right, :)
        temp_point = points(i + 1, :)
        points(i + 1, :) = points(right, :)
        points(right, :) = temp_point

        ! Swap indices(i+1) and indices(right)
        temp_index = indices(i + 1)
        indices(i + 1) = indices(right)
        indices(right) = temp_index

        ! Return the pivot index
        pivot_index = i + 1
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end function partition
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !k-nearest neighbor search
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function knn_search_init(node, target, k) result(query)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        implicit none
        !in
        integer, intent(in) :: k ! Number of nearest neighbors to find
        type(KDTreeNode), pointer, intent(in) :: node
        real, intent(in) :: target(3)
        !local
        integer :: init_depth = 0
        !out
        real :: dist(k)
        integer(kind=8) :: idx(k)
        type(KDTreeResult) :: query

        !Initialize 
        dist = HUGE(0.0)
        idx = -1

        call knn_search(node, init_depth, target, dist, idx, k)

        query%idx = idx
        query%dist = dist

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    endfunction knn_search_init
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    recursive subroutine knn_search(node, depth, target, dist, idx, k)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        implicit none
        !in
        integer :: k ! Number of nearest neighbors to find
        type(KDTreeNode), pointer, intent(in) :: node ! Starting node (usually the root)
        real, intent(in) :: target(3)                 ! Target point (3D)
        integer, intent(in) :: depth     
        !local
        real, intent(inout) :: dist(k)  
        integer(kind=8), intent(inout) :: idx(k) 
        integer :: i
        real :: dist_current, dist_kth
        integer :: axis
        ! Temporary point for contiguous memory access
        real :: temp_point(3)

        if (.not. associated(node)) return

        ! First, check if it is a leaf node
        if (node%is_leaf == 1) then
            
            ! Check all points in the leaf
            do i = 1, size(node%leaf_indices)
                temp_point = node%leaf_points(i, :)
                dist_current = distance(temp_point, target)
                dist_kth = dist(k)

                ! If the current point is closer than the k-th best, update the list
                if (dist_current < dist_kth) then
                    dist(k) = dist_current
                    idx(k) = node%leaf_indices(i)
                    call shift_knn(dist, idx, k)
                end if
            end do

        else

            ! Calculate distances
            dist_current = distance(node%point, target)
            dist_kth = dist(k)

            ! Update best points and indices if the current node is closer than the k-th best
            if (dist_current < dist_kth) then
                dist(k) = dist_current
                idx(k) = node%index
                call shift_knn(dist, idx, k)
            end if

            axis = mod(depth, 3)  ! Determine the current axis (0 for x, 1 for y, 2 for z)

            ! Recursively search the subtree that contains the target
            if (target(axis+1) < node%point(axis+1)) then
                call knn_search(node%left, depth + 1, target, dist, idx, k)
                dist_kth = dist(k)
                !Check if we need to search the right subtree 
                !(dist_kth is still bigger than the distance to the splitting plane)
                if (abs(target(axis+1) - node%point(axis+1)) < 2.*dist_kth) then
                    call knn_search(node%right, depth + 1, target, dist, idx, k)
                end if
            else
                call knn_search(node%right, depth + 1, target, dist, idx, k)
                dist_kth = dist(k)
                !Check if we need to search the right subtree (dist_kth is still bigger than the distance to the splitting plane)
                if (abs(target(axis+1) - node%point(axis+1)) < 2.*dist_kth) then
                    call knn_search(node%left, depth + 1, target, dist, idx, k)
                end if
            end if
        end if

        contains

            !this is only a shiftdown of the k-th element
            !not a full sort
            subroutine shift_knn(dist, idx, k)
                implicit none
                !inout
                real, intent(inout) :: dist(k)
                integer(kind=8), intent(inout) :: idx(k)
                integer, intent(in) :: k
                !local
                real :: temp_dist
                integer(kind=8) :: temp_idx
                integer :: i
            
                ! Store the new element to be inserted
                temp_dist = dist(k)
                temp_idx = idx(k)
            
                ! Start from the end of the array and move the new element to its correct position
                i = k - 1
            
                ! Shift elements greater than temp_dist to the right
                do while (i >= 1)
                    if (dist(i) <= temp_dist) exit  ! Exit the loop if the correct position is found
                    dist(i + 1) = dist(i)
                    idx(i + 1) = idx(i)
                    i = i - 1
                end do
            
                ! Insert the new element into the correct position
                dist(i + 1) = temp_dist
                idx(i + 1) = temp_idx
            end subroutine shift_knn

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine knn_search
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Search for points within a given radius
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function ball_search_init(node, target, radius) result(query)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    implicit none
    !in
    real :: radius ! Radius of the ball
    type(KDTreeNode), pointer, intent(in) :: node
    real, intent(in) :: target(3)
    !local
    integer :: init_depth = 0
    !out
    real, allocatable :: dist(:) ! Distance of the points within the radius
    integer(kind=8), allocatable :: idx(:) !index of the points within the radius
    type(KDTreeResult) :: query

    call ball_search(node, init_depth, target, dist, idx, radius)

    ! No results found
    if (.not. allocated(dist)) then
        allocate(dist(0))
        allocate(idx(0))
    else
        ! Last step, sort the distances
        call quicksort(dist, idx, size(idx))
    end if

    query%idx = idx
    query%dist = dist
    
    contains

        ! Quicksort to sort the distances and indices
        subroutine quicksort(dist, idx, n)
            implicit none
            !in/out
            real, intent(inout) :: dist(n)    ! Distance array to be sorted
            integer(kind=8), intent(inout) :: idx(n) ! Corresponding indices
            integer, intent(in) :: n          ! Number of elements to sort
            !local
            integer :: low, high
        
            low = 1
            high = n
            call quicksort_recursive(dist, idx, low, high, n)
        end subroutine quicksort    
        
        recursive subroutine quicksort_recursive(dist, idx, low, high, n)
            implicit none
            !in/out
            integer, intent(in) :: n
            real, intent(inout) :: dist(n)
            integer(kind=8), intent(inout) :: idx(n)
            integer, intent(in) :: low, high
            !local
            integer :: pivot_index

            if (low < high) then
                ! Partition the array and get the pivot index
                call partition2(dist, idx, low, high, pivot_index, n)

                ! Recursively sort the subarrays
                call quicksort_recursive(dist, idx, low, pivot_index - 1, n)
                call quicksort_recursive(dist, idx, pivot_index + 1, high, n)
            end if
            
        end subroutine quicksort_recursive
        
        subroutine partition2(dist, idx, low, high, pivot_index, n)
            implicit none
            !in/out
            integer, intent(in) :: n
            real, intent(inout) :: dist(n)
            integer(kind=8), intent(inout) :: idx(n)
            integer, intent(in) :: low, high
            integer, intent(out) :: pivot_index
            !local
            real :: pivot_value
            integer :: i, j
            real :: temp_dist
            integer(kind=8) :: temp_idx

            ! Choose the pivot (here, we use the last element)
            pivot_value = dist(high)
            i = low - 1

            ! Partition the array
            do j = low, high - 1
                if (dist(j) <= pivot_value) then
                    i = i + 1
                    ! Swap dist(i) and dist(j)
                    temp_dist = dist(i)
                    dist(i) = dist(j)
                    dist(j) = temp_dist
                    ! Swap idx(i) and idx(j)
                    temp_idx = idx(i)
                    idx(i) = idx(j)
                    idx(j) = temp_idx
                end if
            end do

            ! Place the pivot in its correct position
            i = i + 1
            temp_dist = dist(i)
            dist(i) = dist(high)
            dist(high) = temp_dist
            temp_idx = idx(i)
            idx(i) = idx(high)
            idx(high) = temp_idx

            ! Return the pivot index
            pivot_index = i
        end subroutine partition2

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end function ball_search_init
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    recursive subroutine ball_search(node, depth, target, dist, idx, radius)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    implicit none
    !in
    real :: radius ! Radius of the ball
    type(KDTreeNode), pointer, intent(in) :: node ! Starting node (usually the root)
    real, intent(in) :: target(3)                 ! Target point (3D)
    integer, intent(in) :: depth     
    !out
    integer(kind=8), allocatable, intent(inout) :: idx(:)  ! Index of the points within the radius
    real, allocatable, intent(inout) :: dist(:)
    !local
    integer :: i
    real :: dist_current
    integer :: axis
    ! Temporary point for contiguous memory access
    real :: temp_point(3)

    if (.not. associated(node)) return

        ! First, check if it is a leaf node
        if (node%is_leaf == 1) then
            
            ! Check all points in the leaf
            do i = 1, size(node%leaf_indices)
                temp_point = node%leaf_points(i, :)
                dist_current = distance(temp_point, target)
                if ( dist_current <= radius ) then
                    ! Append the index to the list
                    call int_add_to_list(idx, node%leaf_indices(i))
                    call real_add_to_list(dist, dist_current)
                end if
            end do

        else

            !Calculate this node distance
            dist_current = distance(node%point, target)

            if (dist_current <= radius) then
                call int_add_to_list(idx, node%index)
                call real_add_to_list(dist, dist_current)
            end if

            axis = mod(depth, 3)  ! Determine the current axis (0 for x, 1 for y, 2 for z)

            ! Recursively search the primary subtree
            if (target(axis+1) < node%point(axis+1)) then
                call ball_search(node%left, depth + 1, target, dist, idx, radius)
                ! Check if we need to search the other subtree
                if (abs(target(axis+1) - node%point(axis+1)) <= 2.*radius) then
                    call ball_search(node%right, depth + 1, target, dist, idx, radius)
                end if
            else
                call ball_search(node%right, depth + 1, target, dist, idx, radius)
                ! Check if we need to search the other subtree
                if (abs(target(axis+1) - node%point(axis+1)) <= 2.*radius) then
                    call ball_search(node%left, depth + 1, target, dist, idx, radius)
                end if
            end if

        endif

    contains
        
        !subroutines to append an element to an array
        subroutine int_add_to_list(indices, new_value)
            implicit none
            integer(kind=8), allocatable, intent(inout) :: indices(:)
            integer(kind=8), intent(in) :: new_value
            integer(kind=8), allocatable :: temp(:)
            integer :: n
        
            if (.not. allocated(indices)) then
                allocate(indices(1))
                indices(1) = new_value
            else
                n = size(indices)
                allocate(temp(n + 1))
                temp(1:n) = indices
                temp(n + 1) = new_value
                call move_alloc(temp, indices)  ! Efficient memory transfer
            end if
        end subroutine int_add_to_list

        subroutine real_add_to_list(dist, new_value)
            implicit none
            real, allocatable, intent(inout) :: dist(:)
            real, intent(in) :: new_value
            real, allocatable :: temp(:)
            integer :: n
        
            if (.not. allocated(dist)) then
                allocate(dist(1))
                dist(1) = new_value
            else
                n = size(dist)
                allocate(temp(n + 1))
                temp(1:n) = dist
                temp(n + 1) = new_value
                call move_alloc(temp, dist)  ! Efficient memory transfer
            end if
        end subroutine real_add_to_list

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine ball_search
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Euclidean distance between two 3D points
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function distance(p1, p2) result(dist)
        implicit none
        real, intent(in) :: p1(3), p2(3)
        real :: dist
        dist = sqrt((p1(1) - p2(1))**2 + (p1(2) - p2(2))**2 + (p1(3) - p2(3))**2)
    end function distance
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!#######################################################
end module kdtree_mod
!#######################################################


!demo
program main
      use kdtree_mod
      use omp_lib
!-------------------------------------------------------
!*     From coretran ****************************************
      use variableKind, only: r64
      use m_allocate, only: allocate
      use m_deallocate, only: deallocate
      use m_KdTree, only: KdTree, KdTreeSearch
      use dargdynamicarray_class, only: dArgDynamicArray
!************************************************************
      implicit none

      type(KDTreeNode), pointer :: root
      type(KDTreeResult) :: query
      real, allocatable :: x(:), y(:), z(:)
      real :: ttarget(3) = [50, 50, 50]  ! Target point in 3D
      integer(kind=8) :: n, i, best_index
      integer :: k
      integer :: ncpu
      real, allocatable :: real_dists(:)
      integer :: t1, t2, trate, tmax, counter
      integer(kind=8), allocatable :: indices(:), indices_coretran(:)
      real, allocatable :: dist(:), dist_coretran(:)
      !coretran
      TYPE(KDTREE) :: tree_coretran
      REAL(R64), ALLOCATABLE :: x_coretran(:), y_coretran(:), z_coretran(:)
      !query
      TYPE (KDTREESEARCH) :: SEARCH
      TYPE (DARGDYNAMICARRAY) :: DA

      ! Get the number of threads
      !$OMP PARALLEL
        !$OMP SINGLE
      ncpu = omp_get_num_threads()
        !$OMP END SINGLE
      !$OMP END PARALLEL
      print *, "Number of threads:", ncpu

      ! Number of points (can be larger than INT*4 limit: 2,147,483,647)
      n = 10000000_8

      allocate(real_dists(n))

      ! Allocate and initialize points (example: random points)
      call system_clock(t1,trate,tmax)
      allocate(x(n), y(n), z(n))
      do i = 1, n
            x(i) = rand() * 100  ! Random x-coordinate
            y(i) = rand() * 100  ! Random y-coordinate
            z(i) = rand() * 100  ! Random z-coordinate
      end do
      CALL system_clock(t2,trate,tmax)
      WRITE(*,*) "Time taken to generate random points:", float(t2 - t1)/1e3, "seconds"

      
      ! Build the KD-Tree
      call system_clock(t1,trate,tmax)
      root => build_kdtree_init(x, y, z)
      call system_clock(t1,trate,tmax)
      WRITE(*,*) "Time taken to build KD-Tree:", float(t1 - t2)/1e3, "seconds"
  
      ! Build the KD-Tree coretran
      x_coretran = x
      y_coretran = y
      z_coretran = z
      call system_clock(t1,trate,tmax)
      tree_coretran = KDTREE(x_coretran, y_coretran, z_coretran)
      call system_clock(t1,trate,tmax)
      WRITE(*,*) "CORETRAN Time taken to build KD-Tree:", float(t1 - t2)/1e3, "seconds"
  

      !KNN TEST
      ! Find the nearest neighbor
      k = 64
      allocate(indices(k))
      call system_clock(t1,trate,tmax)
      do i = 1, 1000000
        query = knn_search_init(root, ttarget, k)
      end do
      call system_clock(t2,trate,tmax)
      WRITE(*,*) "Time taken to find nearests neighbors:", float(t2 - t1)/1e3, "seconds"
      indices = query%idx
      dist = query%dist

      !CORETRAN KNN TEST
      call system_clock(t1,trate,tmax)
      do i = 1, 1000000
        DA = SEARCH%KNEAREST(tree_coretran, x_coretran, y_coretran, z_coretran, &
                XQUERY=DBLE(ttarget(1)), YQUERY=DBLE(ttarget(2)),  &
                ZQUERY=DBLE(ttarget(3)), K=k)
      end do
      call system_clock(t2,trate,tmax)
      indices_coretran = DA%I%VALUES
      dist_coretran = DA%V%VALUES
      WRITE(*,*) "CORETRAN Time taken to find nearests neighbors:", float(t2 - t1)/1e3, "seconds"


    !   !Check if all values are in 
    !   ! Print distances to see they are sorted
    !   do i=1,k 
    !     write(*,*) indices(i), dist(i), indices_coretran(i), dist_coretran(i)
    !   end do

    ! ! BALL SEARCH TEST
    !   real_dists = sqrt((x - ttarget(1))**2 + (y - ttarget(2))**2 + (z - ttarget(3))**2)
    !   call system_clock(t1,trate,tmax)
    !   do i=1,1000000
    !     query = ball_search_init(root, ttarget, 0.5)
    !   end do
    !   call system_clock(t2,trate,tmax)
    !   WRITE(*,*) "Time taken to find points within the ball:", float(t2 - t1)/1e3, "seconds"
    !   indices = query%idx
    !   dist = query%dist
    !   print *, "Points within the ball:", size(indices, 1)

    !   call system_clock(t1,trate,tmax)
    !   do i=1,1000000
    !   DA = SEARCH%KNEAREST(tree_coretran, x_coretran, y_coretran, z_coretran, &
    !                XQUERY=DBLE(ttarget(1)), YQUERY=DBLE(ttarget(2)),  &
    !                ZQUERY=DBLE(ttarget(3)), RADIUS=DBLE(0.5))  
    !   end do
    !   call system_clock(t2,trate,tmax)
    !   WRITE(*,*) "CORETRAN Time taken to find points within the ball:", float(t2 - t1)/1e3, "seconds"
    !   indices_coretran = DA%I%VALUES
    !   dist_coretran = DA%V%VALUES
    !   print *, "Coretran points within the ball:", size(indices_coretran)

    !   do i=1,size(indices)
    !     write(*,*) indices(i), dist(i), indices_coretran(i), dist_coretran(i)
    !   end do


      ! Deallocate memory
      deallocate(indices)
      deallocate(x, y, z)
      deallocate(real_dists)
      deallocate(root)

  end program main
