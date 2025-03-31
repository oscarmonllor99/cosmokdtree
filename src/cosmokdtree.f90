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
    public :: build_kdtree, KDTreeNode, KDTreeResult, knn_search, ball_search
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

#if periodic == 1
    ! Box size (periodic boundary conditions), global variable
    real :: L(3) ! Will be initialized in build_kdtree
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
        real, intent(in) :: x(:), y(:), z(:)
#if periodic == 1
        real, intent(in) :: Lx, Ly, Lz
        integer :: flag_stop
#endif
        real, allocatable :: points(:, :)
        integer(kind=8), allocatable :: indices(:)
        integer(kind=8) :: n, i
        integer :: depth, max_depth, nproc, leafsize
        type(KDTreeNode), pointer :: tree
    
        ! Enable nested parallelism
        call omp_set_nested(.true.) 
    
        ! Number of points
        n = size(x, kind=8)

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
    !$OMP PARALLEL
    !$OMP SINGLE

    !$OMP TASK
    node%left => build_kdtree_recursive(points(1:median-1,:),indices(1:median-1),depth+1,max_depth,leafsize)
    !$OMP END TASK

    !$OMP TASK
    node%right => build_kdtree_recursive(points(median+1:,:),indices(median+1:),depth+1,max_depth,leafsize)
    !$OMP END TASK

    !$OMP END SINGLE
    !$OMP END PARALLEL
    else
    node%left => build_kdtree_recursive(points(1:median-1,:),indices(1:median-1),depth+1,max_depth,leafsize)
    node%right => build_kdtree_recursive(points(median+1:,:),indices(median+1:),depth+1,max_depth,leafsize)
    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end function build_kdtree_recursive
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Quickselect function to find the k-th smallest element along a specified axis
    ! Ensures all points below k are less than or equal to the k-th point
    ! and all points above k are greater than or equal to the k-th point with the
    ! specified axis.
    subroutine quickselect(points, indices, k, axis, kth_point, kth_index)
        implicit none
        real, intent(inout) :: points(:, :)      ! 2D array of points
        integer(kind=8), intent(inout) :: indices(:) ! 1D array of indices
        integer(kind=8), intent(in) :: k         ! k-th smallest element to find
        integer, intent(in) :: axis             ! Axis to sort along (0 for x, 1 for y, 2 for z)
        real, intent(out) :: kth_point(size(points, 2))  ! The k-th smallest point
        integer(kind=8), intent(out) :: kth_index        ! Index of the k-th smallest point
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
    end subroutine quickselect

    ! Partition function for Quickselect
    function partition(points, indices, left, right, axis) result(pivot_index)
        implicit none
        real, intent(inout) :: points(:, :)      ! 2D array of points
        integer(kind=8), intent(inout) :: indices(:) ! 1D array of indices
        integer(kind=8), intent(in) :: left, right  ! Left and right bounds of the partition
        integer, intent(in) :: axis             ! Axis to sort along (0 for x, 1 for y, 2 for z)
        integer(kind=8) :: pivot_index, i, j
        real :: pivot_value, temp_point(size(points, 2))
        integer(kind=8) :: temp_index
    
        ! Choose pivot as the middle element
        pivot_index = (left + right) / 2
        pivot_value = points(pivot_index, axis + 1)
    
        ! Move pivot to the end
        call swap(points, indices, pivot_index, right)
    
        i = left - 1
    
        ! Partition the array
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
        real, intent(inout) :: points(:, :)
        integer(kind=8), intent(inout) :: indices(:)
        integer(kind=8), intent(in) :: i, j
        real :: temp_point(size(points, 2))
        integer(kind=8) :: temp_index
    
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
        real, intent(in) :: targett(3)
        !local
        integer :: init_depth = 0
        !out
        real :: dist(k)
        integer(kind=8) :: idx(k)
        type(KDTreeResult) :: query

        !Initialize 
        dist = HUGE(0.0)
        idx = -1

        call knn_search_recursive(node, init_depth, targett, dist, idx, k)

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
        real, intent(in) :: targett(3)                 ! Target point (3D)
        integer, intent(in) :: depth     
        !local
        real, intent(inout) :: dist(k)  
        integer(kind=8), intent(inout) :: idx(k) 
        integer :: i
        real :: dist_current, dist_kth, d1d
        real :: epsilon = 1.e-6
        integer :: axis
        integer :: look_opposite
        ! Temporary point for contiguous memory access
        real :: temp_point(3)

        if (.not. associated(node)) return

        ! First, check if it is a leaf node
        if (node%is_leaf == 1) then
            
            ! Check all points in the leaf
            do i = 1, size(node%leaf_indices)
                temp_point = node%leaf_points(i, :)
                dist_current = distance(temp_point, targett)
                dist_kth = dist(k)

                ! If the current point is closer than the k-th best, update the list
                if (dist_current < dist_kth + epsilon) then
                    dist(k) = dist_current
                    idx(k) = node%leaf_indices(i)
                    call shift_knn(dist, idx, k)
                end if
            end do

        else

            ! Calculate distances
            dist_current = distance(node%point, targett)
            dist_kth = dist(k)

            ! Update best points and indices if the current node is closer than the k-th best
            if (dist_current < dist_kth + epsilon) then
                dist(k) = dist_current
                idx(k) = node%index
                call shift_knn(dist, idx, k)
            end if

            axis = node%axis
            ! 1D distance from target to the splitting plane
            d1d = targett(axis+1) - node%point(axis+1)
            ! Recursively search the subtree that contains the target
            if (d1d < 0) then
                call knn_search_recursive(node%left, depth + 1, targett, dist, idx, k)
                dist_kth = dist(k)
                !Check if we need to search the right subtree 
                look_opposite = 0
                if (abs(d1d) < dist_kth) look_opposite = 1
#if periodic == 1
                if (targett(axis+1) - dist_kth <= -L(axis+1) / 2. ) look_opposite = 1
#endif
                if (look_opposite == 1) then
                    call knn_search_recursive(node%right, depth + 1, targett, dist, idx, k)
                end if
            else
                call knn_search_recursive(node%right, depth + 1, targett, dist, idx, k)
                dist_kth = dist(k)
                !Check if we need to search the left subtree
                if (abs(d1d) < dist_kth) look_opposite = 1
#if periodic == 1
                if (targett(axis+1) + dist_kth >= L(axis+1) / 2. ) look_opposite = 1
#endif
                if (look_opposite == 1) then
                    call knn_search_recursive(node%left, depth + 1, targett, dist, idx, k)
                end if
            end if


        end if !node%is_leaf

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
    end subroutine knn_search_recursive
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Search for points within a given radius
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function ball_search(node, targett, radius) result(query)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    implicit none
    !in
    real :: radius ! Radius of the ball
    type(KDTreeNode), pointer, intent(in) :: node
    real, intent(in) :: targett(3)
    !local
    integer :: init_depth = 0
    integer :: count_idx, count_dist, count ! Counters for the number of elements in idx and dist
    integer(kind=8), allocatable :: temp_idx(:)
    real, allocatable :: temp_dist(:)
    !out
    real, allocatable :: dist(:) ! Distance of the points within the radius
    integer(kind=8), allocatable :: idx(:) !index of the points within the radius
    type(KDTreeResult) :: query

    !Preallocate dist and idx
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
        ! Last step, sort the distances
        call quicksort(dist, idx, size(idx))
    end if

    query%idx = idx
    query%dist = dist

    deallocate(dist, idx)
    
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
    end function ball_search
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
    recursive subroutine ball_search_recursive(node, depth, targett, dist, idx, radius, count_idx, count_dist)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
        implicit none
        !in
        real :: radius ! Radius of the ball
        type(KDTreeNode), pointer, intent(in) :: node ! Starting node (usually the root)
        real, intent(in) :: targett(3)                 ! Target point (3D)
        integer, intent(in) :: depth     
        !out
        integer(kind=8), allocatable, intent(inout) :: idx(:)  ! Index of the points within the radius
        real, allocatable, intent(inout) :: dist(:)
        integer, intent(inout) :: count_idx, count_dist ! Counters for the number of elements in idx and dist
        !local
        integer :: i
        real :: dist_current, d1d
        integer :: axis
        real :: epsilon = 1.e-6
        integer :: look_opposite
        ! Temporary point for contiguous memory access
        real :: temp_point(3)

        if (.not. associated(node)) then
            return
        end if

        ! First, check if it is a leaf node
        if (node%is_leaf == 1) then
            
            ! Check all points in the leaf
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
                look_opposite = 0
                if (abs(d1d) <= radius) look_opposite = 1
#if periodic == 1
                if (targett(axis+1) - radius <= -L(axis+1) / 2. ) look_opposite = 1
#endif
                if (look_opposite == 1) then
                    call ball_search_recursive(node%right, depth + 1, targett, dist, idx, radius, count_idx, count_dist)
                end if
            else
                call ball_search_recursive(node%right, depth + 1, targett, dist, idx, radius, count_idx, count_dist)
                look_opposite = 0
                if (abs(d1d) <= radius) look_opposite = 1
#if periodic == 1
                if (targett(axis+1) + radius >= L(axis+1) / 2. ) look_opposite = 1
#endif
                if (look_opposite == 1) then
                    call ball_search_recursive(node%left, depth + 1, targett, dist, idx, radius, count_idx, count_dist)
                end if
            end if

        endif

    contains

        !subroutines to append an element to an array
        subroutine int_add_to_list(indices, new_value, count)
            implicit none
            integer(kind=8), allocatable, intent(inout) :: indices(:)
            integer(kind=8), intent(in) :: new_value
            integer, intent(inout) :: count
            integer(kind=8), allocatable :: temp(:)
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
                    allocate(temp(10 * n))
                    temp(1:n) = indices
                    call move_alloc(temp, indices)
                end if
                count = count + 1
                indices(count) = new_value
            end if
        end subroutine int_add_to_list

        subroutine real_add_to_list(dist, new_value, count)
            implicit none
            real, allocatable, intent(inout) :: dist(:)
            real, intent(in) :: new_value
            integer, intent(inout) :: count
            real, allocatable :: temp(:)
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
                    allocate(temp(10 * n))
                    temp(1:n) = dist
                    call move_alloc(temp, dist)
                end if
                count = count + 1
                dist(count) = new_value
            end if
        end subroutine real_add_to_list

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine ball_search_recursive
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Euclidean distance between two 3D points
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function distance(p1, p2) result(dist)
        implicit none
        real, intent(in) :: p1(3), p2(3)
        real :: dist
        real :: dx, dy, dz
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

!#######################################################
end module cosmokdtree
!#######################################################