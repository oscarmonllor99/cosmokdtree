program main
    use cosmokdtree
    use omp_lib
!-------------------------------------------------------

    implicit none

    type(KDTreeNode), pointer :: root
    type(KDTreeResult) :: query
    real, allocatable :: x(:), y(:), z(:)
    real :: ttarget(3) = [2., -49.99, 10.]  ! Target point in 3D
    integer(kind=8) :: n, i
    integer :: k
    integer :: ncpu
    integer :: t1, t2, trate, tmax, counter
    integer(kind=8), allocatable :: indices(:)
    real, allocatable :: dist(:)
    real, allocatable :: real_dists(:), temp(:)
    integer(kind=8), allocatable :: real_index(:)
    real :: dx,dy,dz,Lx,Ly,Lz,radius

    ! Get the number of threads
    !$OMP PARALLEL
      !$OMP SINGLE
    ncpu = omp_get_num_threads()
      !$OMP END SINGLE
    !$OMP END PARALLEL
    print *, "Number of threads:", ncpu

    ! Number of points (can be larger than INT*4 limit: 2,147,483,647)
    n = 10000000_8

    allocate(real_dists(n), temp(n))

    ! Allocate and initialize points (example: random points)
    call system_clock(t1,trate,tmax)
    allocate(x(n), y(n), z(n))

    Lx = 100.
    Ly = 100.
    Lz = 100.
    do i = 1, n
          x(i) = (rand()-0.5) * Lx  ! Random x-coordinate
          y(i) = (rand()-0.5) * Ly  ! Random y-coordinate
          z(i) = (rand()-0.5) * Lz  ! Random z-coordinate
    end do
    CALL system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to generate random points:", float(t2 - t1)/1e3, "seconds"
    
    ! Build the KD-Tree
    call system_clock(t1,trate,tmax)
    root => build_kdtree(x, y, z, 100., 100., 100.)
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to build KD-Tree:", float(t2 - t1)/1e3, "seconds"

    real_dists = 0
    do i = 1, n
      dx = x(i) - ttarget(1)
      dy = y(i) - ttarget(2)
      dz = z(i) - ttarget(3)
      dx = min(abs(dx), Lx - abs(dx))
      dy = min(abs(dy), Ly - abs(dy))
      dz = min(abs(dz), Lz - abs(dz))
      real_dists(i) = sqrt(dx**2 + dy**2 + dz**2)
    end do

    !KNN TEST
    ! Find the nearest neighbor
    k = 1000
    allocate(indices(k))
    call system_clock(t1,trate,tmax)
    query = knn_search(root, ttarget, k)
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to find nearest neighbours:", float(t2 - t1)/1e3, "seconds"
    indices = query%idx
    dist = query%dist

    ! !CALCULATING EXACT NEAREST NEIGHBOURS
    ! allocate(real_index(k))
    ! temp = real_dists
    ! do i = 1, k
    !   real_index(i) = minloc(temp, 1)
    !   temp(real_index(i)) = HUGE(0.)
    ! end do

    ! !COMPARING THE RESULTS
    ! do i=1,k
    !   WRITE(*,*) indices(i), real_index(i), dist(i), real_dists(real_index(i))
    ! end do
    ! deallocate(real_index)

    deallocate(indices,dist)


    !BALL SEARCH TEST
    radius = 10.

    call system_clock(t1,trate,tmax)
    query = ball_search(root, ttarget, radius)
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to find points within the ball:", float(t2 - t1)/1e3, "seconds"
    indices = query%idx
    dist = query%dist

    ! !CALCULATING EXACT BALL SEARCH
    ! counter = 0
    ! do i = 1, n
    !   if (real_dists(i) <= radius) then
    !     counter = counter + 1
    !   end if
    ! end do

    ! allocate(real_index(counter))

    ! temp = real_dists
    ! do i = 1, counter
    !   real_index(i) = minloc(temp, 1)
    !   temp(real_index(i)) = HUGE(0.)
    ! end do

    ! !COMPARING THE RESULTS
    ! do i=1,counter
    !   WRITE(*,*) indices(i), real_index(i), dist(i), real_dists(real_index(i))
    ! end do

    ! Deallocate memory
    if (allocated(indices)) deallocate(indices)
    if (allocated(dist)) deallocate(dist)
    if (allocated(real_dists)) deallocate(real_dists)
    if (allocated(temp)) deallocate(temp)
    if (allocated(real_index)) deallocate(real_index)
    deallocate(x, y, z)
    deallocate(root)

end program main
