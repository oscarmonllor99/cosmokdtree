!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! kd-tree module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! GRUP DE COSMOLOGIA COMPUTACIONAL (GCC) UNIVERSITAT DE VALÈNCIA
! Authors: Óscar Monllor Berbegal and David Vallés Pérez
! Date: 30/01/2025
! Last update: 04/05/2025
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Brief description:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! - This module implements an optimized parallel k-d tree construction and 
!   search for k-nearest neighbors and points within a given radius.
!
! - Quickselect is used to find the point that splits the space
!   in two halves along a given axis (median). Use of median ensures 
!   that the tree is balanced.
!
! - The tree is built recursively, with the splitting axis changing
!   at each level (x, y, z, x, y, z, ...). The tree is built in parallel 
!   using OpenMP tasks due to Divide and Conquer nature of the algorithm.
!
! - The search for k-nearest neighbors uses an insertion shiftdown 
!   too quickly sort and replace the nearest neighbours found.
!
! - quicksort is used to sort the distances and indices of    
!   points within a given radius
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Pending improvements:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! - knn-search relies on shiftdown to sort the distances and indices, which
!   has O(k) complexity. This can be improved by using a priority queue (heap).
!
! - Makefile
!
! - Periodic boundary conditions
!
! - Better parallelism for tree building (if possible)
!
! - Better neighbour search (e.g. binary search, or heap...)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!#######################################################
module cosmokdtree
!#######################################################
    implicit none
    private  
    public :: build_kdtree, KDTreeNode, KDTreeResult, knn_search, ball_search, &
               box_search

    !+++++++++++++++++++++++++++++++
    !++++ Precision and integer kind
    !+++++++++++++++++++++++++++++++
#if doubleprecision == 1 
    integer, parameter :: prec = 8 
#else
    integer, parameter :: prec = 4
#endif

#if longint == 1
    integer, parameter :: intkind = 8
#else
    integer, parameter :: intkind = 4
#endif

    !+++++++++++++++++++++++++++++++
    !++++ Type definitions
    !+++++++++++++++++++++++++++++++
    type :: KDTreeNode
        !basic ---------------------
        real(kind=prec):: point(3)       ! 3D point (x, y, z)
        integer(kind=intkind) :: index       ! Index of the point within the original array
        integer :: axis        ! Splitting axis (0 for x, 1 for y, 2 for z)
        type(KDTreeNode), pointer :: left => null()  ! Left child
        type(KDTreeNode), pointer :: right => null() ! Right child
        !leaf, for faster search and building
        integer :: is_leaf     ! Flag to indicate if the node is a leaf
        real(kind=prec), pointer :: leaf_points(:, :) => null()  ! Points in the leaf (for leaf nodes)
        integer(kind=intkind), pointer :: leaf_indices(:) => null()  ! Indices of points in the leaf
    end type KDTreeNode

    type :: KDTreeResult
        integer(kind=intkind), allocatable :: idx(:)
        real(kind=prec), allocatable :: dist(:)
    end type KDTreeResult
    !+++++++++++++++++++++++++++++++

#if periodic == 1
    ! Box size (periodic boundary conditions), global variable
    real(kind=prec):: L(3) ! Will be initialized in build_kdtree
#endif

contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Initialize kd-tree construction
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if periodic == 1
    function build_kdtree(x, y, z, Lx, Ly, Lz) result(tree)
#else
    function build_kdtree(x, y, z) result(tree)
#endif
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        use omp_lib
        implicit none
        real(kind=prec), intent(in) :: x(:), y(:), z(:)
#if periodic == 1
        real(kind=prec), intent(in) :: Lx, Ly, Lz
        integer :: flag_stop
#endif
        real(kind=prec), allocatable :: points(:, :)
        integer(kind=intkind), allocatable :: indices(:)
        integer(kind=intkind) :: n, i
        integer :: depth, max_depth, nproc, leafsize
        type(KDTreeNode), pointer :: tree
    
        ! Enable nested parallelism
        call omp_set_nested(.true.) 
    
        ! Number of points
        n = size(x, kind=intkind)

# if periodic == 1
        ! Set the box size (global variable of the module)
        L = [Lx, Ly, Lz]

        ! Check all points are within (-Lx/2, Lx/2)
        flag_stop = 0
        !$OMP PARALLEL SHARED(x,y,z,n,Lx,Ly,Lz), &
        !$OMP PRIVATE(i)
        !$OMP DO REDUCTION(+:flag_stop)
        do i=1,n
            if(x(i) .lt. -Lx/2 .or. x(i) .gt. Lx/2) flag_stop = 1
            if(y(i) .lt. -Ly/2 .or. y(i) .gt. Ly/2) flag_stop = 1
            if(z(i) .lt. -Lz/2 .or. z(i) .gt. Lz/2) flag_stop = 1
        enddo
        !$OMP ENDDO
        !$OMP END PARALLEL

        if (flag_stop .gt. 0) STOP 'Points outside (-L/2, L/2) range !!'
# endif

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
        max_depth = compute_max_depth(omp_get_max_threads())
        
        ! Leafsize scaling with the number of points
        leafsize = int(real(n)**0.333 / 4.)
        leafsize = max(leafsize, 1)

        tree => build_kdtree_recursive(points, indices, depth, max_depth, leafsize)

        deallocate(points, indices)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end function build_kdtree
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
    recursive function build_kdtree_recursive(points,indices,depth,max_depth,leafsize) result(node)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    use omp_lib
    implicit none
    real(kind=prec), intent(inout) :: points(:, :)      ! 2D array of points
    integer(kind=intkind), intent(inout) :: indices(:) ! 1D array of indices
    integer, intent(in) :: depth         ! Current depth in the tree
    integer, intent(in) :: max_depth ! max_depth to allow parallelism
    integer, intent(in) :: leafsize ! size of leaf nodes
    integer :: axis !axis to split points
    real(kind=prec) :: var(3) ! variance of the points in each axis
    type(KDTreeNode), pointer :: node   ! New node to be created
    real(kind=prec):: kth_point(size(points, 2)) ! median point
    integer(kind=intkind) :: median, median_approx ! half size of data
    integer(kind=intkind) :: kth_index ! index of the median point

    if (size(points, 1) == 0) then
        node => null()
        return
    end if

    ! Splitting axis: maximum variance
    var(1) = variance(points, 0) ! variance in x
    var(2) = variance(points, 1) ! variance in y
    var(3) = variance(points, 2) ! variance in z
    axis = 0
    if (var(2) > var(1)) axis = 1
    if (var(3) > var(2)) axis = 2

    ! Find median and partition points
    median = size(points, 1, kind=intkind) / 2 + 1
    call quickselect(points, indices, median, axis, kth_point, kth_index, median_approx)

    ! Allocate node
    allocate(node)
    node%point = kth_point
    node%index = kth_index
    node%axis = axis

    ! Check if this is a leaf node
    if (size(points, 1, kind=intkind) <= leafsize) then
        node%is_leaf = 1
        ! Store all points and indices in the leaf node
        allocate(node%leaf_points(size(points, 1, kind=intkind), size(points, 2)))
        allocate(node%leaf_indices(size(points, 1, kind=intkind)))
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
    !$OMP PARALLEL
    !$OMP SINGLE

    !$OMP TASK
    node%left => build_kdtree_recursive(points(1:median_approx-1,:),indices(1:median_approx-1),depth+1,max_depth,leafsize)
    !$OMP END TASK

    !$OMP TASK
    node%right => build_kdtree_recursive(points(median_approx+1:,:),indices(median_approx+1:),depth+1,max_depth,leafsize)
    !$OMP END TASK

    !$OMP END SINGLE
    !$OMP END PARALLEL
    else
    node%left => build_kdtree_recursive(points(1:median_approx-1,:),indices(1:median_approx-1),depth+1,max_depth,leafsize)
    node%right => build_kdtree_recursive(points(median_approx+1:,:),indices(median_approx+1:),depth+1,max_depth,leafsize)
    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end function build_kdtree_recursive
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function variance(points, axis) result(var)
        implicit none
        real(kind=prec), intent(in) :: points(:, :)      ! 2D array of points
        integer, intent(in) :: axis             ! Axis to calculate variance (0 for x, 1 for y, 2 for z)
        real(kind=prec) :: var, sum, sum_sq
        integer(kind=intkind) :: n, i

        sum = 0.
        sum_sq = 0.
        n = size(points, 1, intkind)

        do i = 1, n
            sum = sum + points(i, axis + 1)
            sum_sq = sum_sq + points(i, axis + 1)**2
        end do

        var = sum_sq / n - (sum / n)**2
    end function variance
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Quickselect function to find the k-th smallest element along a specified axis
    ! Ensures all points below k are less than or equal to the k-th point
    ! and all points above k are greater than or equal to the k-th point with the
    ! specified axis.
    subroutine quickselect(points, indices, k, axis, kth_point, kth_index, k_exit)
        implicit none
        !in
        real(kind=prec), intent(inout) :: points(:, :)      ! 2D array of points
        integer(kind=intkind), intent(inout) :: indices(:) ! 1D array of indices
        integer(kind=intkind), intent(in) :: k         ! k-th smallest element to find
        integer, intent(in) :: axis             ! Axis to sort along (0 for x, 1 for y, 2 for z)
        !local
        integer(kind=intkind) :: left, right, pivot_index
        integer(kind=intkind) :: k_up, k_down
        !out
        real(kind=prec), intent(out) :: kth_point(size(points, 2))  ! The k-th smallest point
        integer(kind=intkind), intent(out) :: kth_index        ! Index of the k-th smallest point
        integer(kind=intkind) :: k_exit
        
        left = 1
        right = size(points, 1, intkind)
        k_up = k + int(0.2 * size(points, 1, intkind))
        k_down = k - int(0.2 * size(points, 1, intkind))
        
        do while (left <= right)
            ! Partition the array and get the pivot index
            pivot_index = partition(points, indices, left, right, axis)
            ! Early exit if the pivot index is within some error margin of k
            if (pivot_index >= k_down .and. pivot_index <= k_up) then
                ! Found the k-th smallest element
                kth_point = points(pivot_index, :)
                kth_index = indices(pivot_index)
                k_exit = pivot_index
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
    end subroutine quickselect

    ! Partition function for Quickselect
    function partition(points, indices, left, right, axis) result(pivot_index)
        implicit none
        real(kind=prec), intent(inout) :: points(:, :)      ! 2D array of points
        integer(kind=intkind), intent(inout) :: indices(:) ! 1D array of indices
        integer(kind=intkind), intent(in) :: left, right  ! Left and right bounds of the partition
        integer, intent(in) :: axis             ! Axis to sort along (0 for x, 1 for y, 2 for z)
        integer(kind=intkind) :: pivot_index, i, j
        real(kind=prec):: pivot_value
    
        ! Choose pivot as the middle element
        pivot_index = (left + right) / 2
        pivot_value = points(pivot_index, axis + 1)
    
        ! Move pivot to the end
        call swap(points, indices, pivot_index, right)
    
        i = left - 1
    
        ! Move all elements smaller than or equal to pivot to the left
        do j = left, right - 1
            if (points(j, axis + 1) <= pivot_value) then
                i = i + 1
                call swap(points, indices, i, j)
            end if
        end do
    
        ! Move pivot to its final position
        call swap(points, indices, i + 1, right)
    
        ! Return the pivot index
        pivot_index = i + 1
    end function partition

    subroutine swap(points, indices, i, j)
        implicit none
        real(kind=prec), intent(inout) :: points(:, :)
        integer(kind=intkind), intent(inout) :: indices(:)
        integer(kind=intkind), intent(in) :: i, j
        real(kind=prec):: temp_point(size(points, 2))
        integer(kind=intkind) :: temp_index
    
        ! Swap points
        temp_point = points(i, :)
        points(i, :) = points(j, :)
        points(j, :) = temp_point
    
        ! Swap indices
        temp_index = indices(i)
        indices(i) = indices(j)
        indices(j) = temp_index
    end subroutine swap
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !k-nearest neighbor search
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function knn_search(node, targett, k) result(query)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        implicit none
        !in
        integer, intent(in) :: k ! Number of nearest neighbors to find
        type(KDTreeNode), pointer, intent(in) :: node
        real(kind=prec), intent(in) :: targett(3)
        !local
        integer :: init_depth = 0
        !out
        real(kind=prec):: dist(k)
        integer(kind=intkind) :: idx(k)
        type(KDTreeResult) :: query

        !Initialize 
        dist = HUGE(0.0)
        idx = -1

        call knn_search_recursive(node, init_depth, targett, dist, idx, k)

        !Perform full sort with quicksort
        call quicksort(dist, idx, k)

        query%idx = idx
        query%dist = dist

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    endfunction knn_search
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    recursive subroutine knn_search_recursive(node, depth, targett, dist, idx, k)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        implicit none
        !in
        integer :: k ! Number of nearest neighbors to find
        type(KDTreeNode), pointer, intent(in) :: node ! Starting node (usually the root)
        real(kind=prec), intent(in) :: targett(3)                 ! Target point (3D)
        integer, intent(in) :: depth     
        !local
        real(kind=prec), intent(inout) :: dist(k)  
        integer(kind=intkind), intent(inout) :: idx(k) 
        integer :: i ! Index running over leaf points
        real(kind=prec):: dist_current, dist_furthest, d1d
        real(kind=prec):: epsilon = 1.e-6 
        integer :: axis
        logical :: look_opposite
        ! Temporary point for contiguous memory access
        real(kind=prec):: temp_point(3)

        if (.not. associated(node)) return

        ! First, check if it is a leaf node
        if (node%is_leaf == 1) then
            
            ! Check all points in the leaf with brute force
            do i = 1, size(node%leaf_indices)
                temp_point = node%leaf_points(i, :)
                dist_current = distance(temp_point, targett)
                dist_furthest = dist(1)

                ! If the current point is closer than the furthest, replace it
                if (dist_current < dist_furthest + epsilon) then
                    dist(1) = dist_current
                    idx(1) = node%leaf_indices(i)
                    call max_heap_insert(dist, idx, k)
                end if
            end do

        else

            ! Calculate distances
            dist_current = distance(node%point, targett)
            dist_furthest = dist(1)

            ! Update best points and indices if the current node is closer than the furthest
            if (dist_current < dist_furthest + epsilon) then
                dist(1) = dist_current
                idx(1) = node%index
                call max_heap_insert(dist, idx, k)
            end if

            axis = node%axis
            ! 1D distance from target to the splitting plane
            d1d = targett(axis+1) - node%point(axis+1)
            ! Recursively search the subtree that contains the target
            if (d1d < 0) then
                call knn_search_recursive(node%left, depth + 1, targett, dist, idx, k)
                dist_furthest = dist(1)
                !Check if we need to search the right subtree 
                look_opposite = .false.
                if (abs(d1d) < dist_furthest) look_opposite = .true.
#if periodic == 1
                if (targett(axis+1) - dist_furthest <= -L(axis+1) / 2. ) look_opposite = .true.
#endif
                if (look_opposite .eqv. .true.) then
                    call knn_search_recursive(node%right, depth + 1, targett, dist, idx, k)
                end if
            else
                call knn_search_recursive(node%right, depth + 1, targett, dist, idx, k)
                dist_furthest = dist(1)
                !Check if we need to search the left subtree
                look_opposite = .false.
                if (abs(d1d) < dist_furthest) look_opposite = .true.
#if periodic == 1
                if (targett(axis+1) + dist_furthest >= L(axis+1) / 2. ) look_opposite = .true.
#endif
                if (look_opposite .eqv. .true.) then
                    call knn_search_recursive(node%left, depth + 1, targett, dist, idx, k)
                end if
            end if


        end if !node%is_leaf

        contains

            !Instead of keeping a sorted list of distances and indices, we use a max-heap
            !to keep track of the furthest point in the list of k nearest neighbors
            subroutine max_heap_insert(dist, idx, k)
                implicit none
                ! Inputs
                integer, intent(in) :: k
                real(kind=prec), intent(inout) :: dist(k)
                integer(kind=intkind), intent(inout) :: idx(k)
                ! Local
                integer :: i, largest, left, right
                real(kind=prec) :: tmp_dist
                integer(kind=intkind) :: tmp_idx
            
                ! The new value has replaced the root value (furthest) at index 1
                ! dist(1), idx(1)

                ! Heapify down from root to restore max-heap property (every parent bigger than its children)
                i = 1
                do
                    left = 2 * i
                    right = 2 * i + 1
                    largest = i
            
                    if (left <= k .and. dist(left) > dist(largest)) largest = left
                    if (right <= k .and. dist(right) > dist(largest)) largest = right
            
                    if (largest /= i) then
                        ! Swap i-th element with largest
                        tmp_dist = dist(i)
                        dist(i) = dist(largest)
                        dist(largest) = tmp_dist
            
                        tmp_idx = idx(i)
                        idx(i) = idx(largest)
                        idx(largest) = tmp_idx
            
                        i = largest

                    else !already fulfills max-heap property
                        exit
                    end if
                end do
            end subroutine max_heap_insert            


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine knn_search_recursive
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Search for points within a given radius
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function ball_search(node, targett, radius) result(query)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    implicit none
    !in
    real(kind=prec):: radius ! Radius of the ball
    type(KDTreeNode), pointer, intent(in) :: node
    real(kind=prec), intent(in) :: targett(3)
    !local
    integer :: init_depth = 0
    integer :: count_idx, count_dist ! Counters for the number of elements in idx and dist
    integer(kind=intkind), allocatable :: temp_idx(:)
    real(kind=prec), allocatable :: temp_dist(:)
    !out
    real(kind=prec), allocatable :: dist(:) ! Distance of the points within the radius
    integer(kind=intkind), allocatable :: idx(:) !index of the points within the radius
    type(KDTreeResult) :: query

    !Preallocate dist and idx with a reasonable size
    !If the size is not enough, it will be resized with 10x the current size
    allocate(dist(1000)) 
    allocate(idx(1000))
    dist = HUGE(0.0)
    idx = -1
    count_dist = 0
    count_idx = 0

    call ball_search_recursive(node, init_depth, targett, dist, idx, radius, count_idx, count_dist)

    !Check
    if (count_idx .ne. count_dist) then
        STOP "Error: count_idx and count_dist do not match!"
    end if

    if (.not. allocated(dist)) STOP 'dist and idx arrays are not allocated!'
        
    if (count_idx == 0) then
        ! No points found
        temp_dist = dist(1:0)
        call move_alloc(temp_dist, dist)
        temp_idx = idx(1:0)
        call move_alloc(temp_idx, idx)

    else
        ! Reallocation to the correct size
        temp_dist = dist(1:count_dist)
        call move_alloc(temp_dist, dist)
        temp_idx = idx(1:count_idx)
        call move_alloc(temp_idx, idx)
        ! Last step, sort the distances (slightly decreases performance)
        call quicksort(dist, idx, size(idx))

    end if

    query%idx = idx
    query%dist = dist

    deallocate(dist, idx)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end function ball_search
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
    recursive subroutine ball_search_recursive(node, depth, targett, dist, idx, radius, count_idx, count_dist)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
        implicit none
        !in
        real(kind=prec):: radius ! Radius of the ball
        type(KDTreeNode), pointer, intent(in) :: node ! Starting node (usually the root)
        real(kind=prec), intent(in) :: targett(3)                 ! Target point (3D)
        integer, intent(in) :: depth     
        !out
        integer(kind=intkind), allocatable, intent(inout) :: idx(:)  ! Index of the points within the radius
        real(kind=prec), allocatable, intent(inout) :: dist(:)
        integer, intent(inout) :: count_idx, count_dist ! Counters for the number of elements in idx and dist
        !local
        integer :: i
        real(kind=prec):: dist_current, d1d
        integer :: axis
        real(kind=prec):: epsilon = 1.e-6
        logical :: look_opposite
        ! Temporary point for contiguous memory access
        real(kind=prec):: temp_point(3)

        if (.not. associated(node)) then
            return
        end if

        ! First, check if it is a leaf node
        if (node%is_leaf == 1) then
            
            ! Check all points in the leaf with brute force
            do i = 1, size(node%leaf_indices)
                temp_point = node%leaf_points(i, :)
                dist_current = distance(temp_point, targett)
                if ( dist_current <= radius + epsilon ) then
                    ! Append the index to the list
                    call int_add_to_list(idx, node%leaf_indices(i), count_idx)
                    call real_add_to_list(dist, dist_current, count_dist)
                end if
            end do

        else

            !Calculate this node distance
            dist_current = distance(node%point, targett)
            if (dist_current <= radius + epsilon) then
                call int_add_to_list(idx, node%index, count_idx)
                call real_add_to_list(dist, dist_current, count_dist)
            end if

            axis = node%axis
            ! 1D distance from target to the splitting plane
            d1d = targett(axis+1) - node%point(axis+1)
            ! Recursively search the subtree that contains the target
            if (d1d < 0) then
                call ball_search_recursive(node%left, depth + 1, targett, dist, idx, radius, count_idx, count_dist)
                ! Check if we need to search the other subtree
                look_opposite = .false.
                if (abs(d1d) <= radius) look_opposite = .true.
#if periodic == 1
                if (targett(axis+1) - radius <= -L(axis+1) / 2. ) look_opposite = .true.
#endif
                if (look_opposite .eqv. .true.) then
                    call ball_search_recursive(node%right, depth + 1, targett, dist, idx, radius, count_idx, count_dist)
                end if
            else
                call ball_search_recursive(node%right, depth + 1, targett, dist, idx, radius, count_idx, count_dist)
                look_opposite = .false.
                if (abs(d1d) <= radius) look_opposite = .true.
#if periodic == 1
                if (targett(axis+1) + radius >= L(axis+1) / 2. ) look_opposite = .true.
#endif
                if (look_opposite .eqv. .true.) then
                    call ball_search_recursive(node%left, depth + 1, targett, dist, idx, radius, count_idx, count_dist)
                end if
            end if

        endif

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine ball_search_recursive
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Search for points within a given box (cuboid) aligned with cartesian axes
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function box_search(node, box) result(query)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        implicit none
        !in
        real(kind=prec) :: box(6) ! Box limits (xmin, xmax, ymin, ymax, zmin, zmax)
        type(KDTreeNode), pointer, intent(in) :: node
        !local
        integer :: init_depth = 0
        integer :: count_idx ! Counters for the number of elements in idx and dist
        integer(kind=intkind), allocatable :: temp_idx(:)
        !out
        integer(kind=intkind), allocatable :: idx(:) !index of the points within the radius
        type(KDTreeResult) :: query
    
        !Preallocate idx with a reasonable size
        !If the size is not enough, it will be resized with 10x the current size
        allocate(idx(1000))
        idx = -1
        count_idx = 0
    
        call box_search_recursive(node, init_depth, idx, box, count_idx)
    
        if (.not. allocated(idx)) STOP 'idx array is not allocated!'
            
        if (count_idx == 0) then
            ! No points found
            temp_idx = idx(1:0)
            call move_alloc(temp_idx, idx)
        else
            ! Reallocation to the correct size
            temp_idx = idx(1:count_idx)
            call move_alloc(temp_idx, idx)
            
            !For the moment, we do not sort the points by distance to the target

        end if
    
        query%idx = idx

        ! Distances are not calculated in this function, hence dist is not allocated

        deallocate(idx)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end function box_search
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
    recursive subroutine box_search_recursive(node, depth, idx, box, count_idx)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
        implicit none
        !in
        real(kind=prec) :: box(6) ! Box limits (xmin, xmax, ymin, ymax, zmin, zmax)
        type(KDTreeNode), pointer, intent(in) :: node ! Starting node (usually the root)
        integer, intent(in) :: depth     
        !out
        integer(kind=intkind), allocatable, intent(inout) :: idx(:)  ! Index of the points within the box
        integer, intent(inout) :: count_idx ! Counters for the number of elements in idx 
        !local
        integer :: i
        logical :: in_box ! tells if a point is inside the query box
        integer :: axis
        ! Temporary point for contiguous memory access
        real(kind=prec):: temp_point(3)

        if (.not. associated(node)) then
            return
        end if

        ! First, check if it is a leaf node
        if (node%is_leaf == 1) then
            
            ! Check all points in the leaf with brute force
            do i = 1, size(node%leaf_indices)
                temp_point = node%leaf_points(i, :)
                call is_in_box(temp_point, box, in_box)
                if ( in_box .eqv. .true. ) then
                    ! Append the index to the list
                    call int_add_to_list(idx, node%leaf_indices(i), count_idx)
                end if
            end do

        else

            !See if this node is inside the box
            call is_in_box(node%point, box, in_box)
            if ( in_box .eqv. .true. ) then
                call int_add_to_list(idx, node%index, count_idx)
            end if

            axis = node%axis

            ! Recursively search the subtrees intersecting the box
            ! Periodicity does not makes sense as the query box is assumed to be contained
            ! inside the (periodic) bounding box of all points
            if (box(2*axis+1) < node%point(axis+1)) then
                call box_search_recursive(node%left, depth + 1, idx, box, count_idx)
            endif

            if (box(2*axis+2) > node%point(axis+1)) then
                call box_search_recursive(node%right, depth + 1, idx, box, count_idx)
            endif

        endif !node%is_leaf

    contains

        !subroutine to check if a point is inside the query box
        subroutine is_in_box(point, box, in_box)
            implicit none
            real(kind=prec), intent(in) :: point(3) ! Point to check
            real(kind=prec), intent(in) :: box(6)   ! Box limits (xmin, xmax, ymin, ymax, zmin, zmax)
            logical, intent(out) :: in_box          ! Result: true if the point is inside the box

            in_box = .true.
            if (point(1) < box(1) .or. point(1) > box(2)) in_box = .false.
            if (point(2) < box(3) .or. point(2) > box(4)) in_box = .false.
            if (point(3) < box(5) .or. point(3) > box(6)) in_box = .false.
        end subroutine is_in_box

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine box_search_recursive
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !subroutines to append an element to an array, called by BALL_SEARCH and BOX_SEARCH
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine int_add_to_list(indices, new_value, count)
        implicit none
        integer(kind=intkind), allocatable, intent(inout) :: indices(:)
        integer(kind=intkind), intent(in) :: new_value
        integer, intent(inout) :: count 
        integer(kind=intkind), allocatable :: temp(:)
        integer :: n

        if (.not. allocated(indices)) then
            ! Initial allocation with a reasonable size
            allocate(indices(1000))
            indices(1) = new_value
            count = 1
        else
            if (count == size(indices)) then
                ! Resize the array
                n = size(indices)
                allocate(temp(2 * n))
                temp(1:n) = indices
                call move_alloc(temp, indices)
            end if
            count = count + 1
            indices(count) = new_value
        end if
    end subroutine int_add_to_list


    subroutine real_add_to_list(dist, new_value, count)
        implicit none
        real(kind=prec), allocatable, intent(inout) :: dist(:)
        real(kind=prec), intent(in) :: new_value
        integer, intent(inout) :: count
        real(kind=prec), allocatable :: temp(:)
        integer :: n

        if (.not. allocated(dist)) then
            ! Initial allocation with a reasonable size
            allocate(dist(1000))
            dist(1) = new_value
            count = 1
        else
            if (count == size(dist)) then
                ! Resize the array by doubling its size
                n = size(dist)
                allocate(temp(2 * n))
                temp(1:n) = dist
                call move_alloc(temp, dist)
            end if
            count = count + 1
            dist(count) = new_value
        end if
    end subroutine real_add_to_list
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Euclidean distance between two 3D points
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function distance(p1, p2) result(dist)
        implicit none
        real(kind=prec), intent(in) :: p1(3), p2(3)
        real(kind=prec):: dist
        real(kind=prec):: dx, dy, dz
        dx = p1(1) - p2(1)
        dy = p1(2) - p2(2)
        dz = p1(3) - p2(3)
#if periodic == 1
        dx = min(abs(dx), L(1) - abs(dx))
        dy = min(abs(dy), L(2) - abs(dy))
        dz = min(abs(dz), L(3) - abs(dz))
#endif
        dist = sqrt(dx**2 + dy**2 + dz**2)
    end function distance
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !Quicksort algorithm for sorting final lists in queries
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine quicksort(dist, idx, n)
        implicit none
        !in/out
        real(kind=prec), intent(inout) :: dist(n)    ! Distance array to be sorted
        integer(kind=intkind), intent(inout) :: idx(n) ! Corresponding indices
        integer, intent(in) :: n          ! Number of elements to sort
        !local
        integer :: low, high
    
        low = 1
        high = n
        call quicksort_recursive(dist, idx, low, high)
    end subroutine quicksort    
    
    recursive subroutine quicksort_recursive(dist, idx, low, high)
        implicit none
        !in/out
        real(kind=prec), intent(inout) :: dist(:)
        integer(kind=intkind), intent(inout) :: idx(:)
        integer, intent(in) :: low, high
        !local
        integer :: pivot_index

        if (low < high) then
            ! Partition the array and get the pivot index
            call partition2(dist, idx, low, high, pivot_index)

            ! Recursively sort the subarrays
            call quicksort_recursive(dist, idx, low, pivot_index - 1)
            call quicksort_recursive(dist, idx, pivot_index + 1, high)
        end if
 
    end subroutine quicksort_recursive

    subroutine partition2(dist, idx, low, high, pivot_index)
        implicit none
        !in/out
        real(kind=prec), intent(inout) :: dist(:)
        integer(kind=intkind), intent(inout) :: idx(:)
        integer, intent(in) :: low, high
        integer, intent(out) :: pivot_index
        !local
        integer :: mid
        real(kind=prec):: pivot_value
        integer :: i, j

        ! Choose the pivot (median of low, mid, high)
        mid = (low + high) / 2

        ! Find median of low, mid, high
        if (dist(low) > dist(mid)) call swap2(dist, idx, low, mid)
        if (dist(low) > dist(high)) call swap2(dist, idx, low, high)
        if (dist(mid) > dist(high)) call swap2(dist, idx, mid, high)
        ! Now dist(low) <= dist(mid) <= dist(high)
        pivot_value = dist(mid)
        ! Move pivot to the end
        call swap2(dist, idx, mid, high)

        i = low - 1

        ! Partition the array
        do j = low, high - 1
            if (dist(j) <= pivot_value) then
                i = i + 1
                call swap2(dist, idx, i, j)
            end if
        end do

        ! Place the pivot in its correct position
        i = i + 1
        call swap2(dist, idx, i, high)

        ! Return the pivot index
        pivot_index = i
    end subroutine partition2

    subroutine swap2(dist, idx, i, j)
        implicit none
        real(kind=prec), intent(inout) :: dist(:)
        integer(kind=intkind), intent(inout) :: idx(:)
        integer(kind=intkind), intent(in) :: i, j
        real(kind=prec):: temp_dist
        integer(kind=intkind) :: temp_idx

        ! Swap dist
        temp_dist = dist(i)
        dist(i) = dist(j)
        dist(j) = temp_dist

        ! Swap idx
        temp_idx = idx(i)
        idx(i) = idx(j)
        idx(j) = temp_idx
    end subroutine swap2
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!#######################################################
end module cosmokdtree
!#######################################################