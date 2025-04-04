program main
    use cosmokdtree
    use omp_lib
!-------------------------------------------------------

    implicit none

    integer, parameter :: prec = 4 
    integer, parameter :: intkind = 4

    type(KDTreeNode), pointer :: root
    type(KDTreeResult) :: query
    real(kind=prec), allocatable :: x(:), y(:), z(:)
    integer(kind=intkind) :: n, i
    integer :: count, count2
    integer :: ncpu
    integer :: t1, t2, trate, tmax
    integer(kind=intkind), allocatable :: indices(:), indices_BF(:)
    real(kind=prec) :: box(6)

    ! Get the number of threads
    !$OMP PARALLEL
      !$OMP SINGLE
    ncpu = omp_get_num_threads()
      !$OMP END SINGLE
    !$OMP END PARALLEL
    print *, "Number of threads:", ncpu

    ! Number of points (can be larger than INT*4 limit: 2,147,483,647)
    n = 50000000
    print *, "Number of points:", n

    ! Allocate and initialize points (example: random points)
    call system_clock(t1,trate,tmax)
    allocate(x(n), y(n), z(n))
    do i = 1, n
          x(i) = (rand()-0.5) * 100  ! Random x-coordinate
          y(i) = (rand()-0.5) * 100  ! Random y-coordinate
          z(i) = (rand()-0.5) * 100  ! Random z-coordinate
    end do
    CALL system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to generate random points:", float(t2 - t1)/1e3, "seconds"

    
    ! Build the KD-Tree
    call system_clock(t1,trate,tmax)
    root => build_kdtree(x, y, z)
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to build KD-Tree:", float(t2 - t1)/1e3, "seconds"

    ! Box query the tree
    ! Query box
    box = [-20., 20., -30., 30., -10., 10.]
    call system_clock(t1,trate,tmax)
    query = box_search(root, box)
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken for box search:", float(t2 - t1)/1e3, "seconds"
    indices = query%idx
    count = size(indices)

    ! Brute-force
    allocate(indices_BF(n))
    count2 = 0
    call system_clock(t1,trate,tmax)
    do i=1,n
      if (x(i) >= box(1) .and. x(i) <= box(2) .and. &
          y(i) >= box(3) .and. y(i) <= box(4) .and. &
          z(i) >= box(5) .and. z(i) <= box(6)) then
        count2 = count2 + 1
        indices_BF(count2) = i
      end if
    end do
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken for brute-force search:", float(t2 - t1)/1e3, "seconds"

    write(*,*) "Number of points in box:", count
    ! write(*,*) "Indices of points in box (KD-Tree):"
    ! write(*,*) indices(1:count)
    write(*,*) "Number of points in box (brute-force):", count2
    ! write(*,*) "Indices of points in box (brute-force):"
    ! write(*,*) indices_BF(1:count2)

    ! Deallocate memory
    deallocate(indices, indices_BF)
    deallocate(x, y, z)
    deallocate(root)

end program main
