program main
!-------------------------------------------------------
!*     From coretran ****************************************
    use variableKind, only: r64
    use m_allocate, only: allocate
    use m_deallocate, only: deallocate
    use m_KdTree, only: KdTree, KdTreeSearch
    use dargdynamicarray_class, only: dArgDynamicArray
!************************************************************

    implicit none

    !parameters
    integer :: ncpu
    integer :: ndim, j
    integer, parameter :: prec = 4 
    integer, parameter :: intkind = 4

    !input data
    real(kind=prec), allocatable :: L(:)
    real(kind=prec), allocatable :: x(:), y(:), z(:)
    real(kind=prec), allocatable :: points(:,:)
    integer(kind=intkind) :: n, i, k
    integer(kind=intkind), allocatable :: ninputs(:)
    

    !time
    integer*8 :: t1, t2, trate, tmax
    real*4, allocatable :: times(:), means(:), stds(:)
    real*4 :: mean, std
    integer ::  nbuild

    !queries
    real(kind=prec), allocatable :: ttarget(:,:)
    real(kind=prec) :: ball_radius

    !coretran
    TYPE(KDTREE) :: tree_coretran
    REAL(R64), ALLOCATABLE :: x_coretran(:), y_coretran(:), z_coretran(:)

    ! Input data
    n = 10000000
    print *, "Number of points:", n, " (", n/1000000, "M)"

    ! Allocate and initialize points (example: random points)
    ! Set the number of dimensions
    ndim = 3
    nbuild = 10
    allocate(L(ndim))
    allocate(ninputs(4))
    allocate(means(4))
    allocate(stds(4))
    ninputs(1) = 1000000
    ninputs(2) = 10000000
    ninputs(3) = 100000000
    ninputs(4) = 1000000000
    means = 0.
    stds = 0.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    do k = 1, size(ninputs)
      n = ninputs(k)
      print *, "Building KD-Tree with ", n, " points."
    
      ! random uniform points in a cube of size L
      allocate(points(n,ndim))
      L = 100.0
      do j = 1, ndim
        do i = 1, n
            points(i,j) = (rand()-0.5) * L(j)  
        end do
      end do

      ! Build the KD-Tree and time it
      allocate( x_coretran(n) )
      allocate( y_coretran(n) )
      allocate( z_coretran(n) )
      x_coretran = points(:,1)
      y_coretran = points(:,2)
      z_coretran = points(:,3)
      allocate(times(nbuild))  
      do i = 1, nbuild
        write(*,*) " Build iteration ", i, " of ", nbuild
        call system_clock(t1,trate,tmax)
        tree_coretran = KDTREE(x_coretran, y_coretran, z_coretran)
        call system_clock(t2,trate,tmax)
        times(i) = float(t2 - t1)/float(trate)
        call tree_coretran%deallocate()
      end do
      
      mean = sum(times)/real(nbuild)
      std = sqrt(sum((times - mean)**2)/real(nbuild))

      means(k) = mean
      stds(k) = std

      deallocate(times)
      deallocate(points)
      deallocate(x_coretran, y_coretran, z_coretran)
    end do
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Save results to a file
    open(unit=10, file='../../../../benchmark/coretran_build_times.txt', status='replace')
    write(10,*) 'N_points  Mean_time(s)  Std_dev(s)'
    do k = 1, size(ninputs)
      write(10,*) ninputs(k), means(k), stds(k)
    end do
    close(10)
    deallocate(ninputs)
    deallocate(L)
    deallocate(means)
    deallocate(stds)  

end program main
