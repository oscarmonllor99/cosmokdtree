program main
    use cosmokdtree
    use omp_lib
!-------------------------------------------------------

    implicit none

    !parameters
    integer :: ncpu
    integer :: ndim, j
    integer, parameter :: prec = 4 
    integer, parameter :: intkind = 4

    !input data
    real(kind=prec), allocatable :: L(:)
    real(kind=prec), allocatable :: points(:,:)
    integer(kind=intkind) :: n, i, k
    integer(kind=intkind), allocatable :: ninputs(:)
    
    !tree results
    type(KDTreeNode), pointer :: root

    !time
    integer*8 :: t1, t2, trate, tmax
    real*4, allocatable :: times(:), means(:), stds(:)
    real*4 :: mean, std
    integer ::  nbuild


    ! Get the number of threads
    !$OMP PARALLEL
      !$OMP SINGLE
    ncpu = omp_get_num_threads()
      !$OMP END SINGLE
    !$OMP END PARALLEL
    print *, "Number of threads:", ncpu

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
    means = 0.0
    stds = 0.0
    
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
      allocate(times(nbuild))  
      do i = 1, nbuild
        write(*,*) " Build iteration ", i, " of ", nbuild
        call system_clock(t1,trate,tmax)
        root => build_kdtree(points)
        call system_clock(t2,trate,tmax)
        times(i) = float(t2 - t1)/float(trate)
        call deallocate_kdtree(root)
      end do
      
      mean = sum(times)/real(nbuild)
      std = sqrt(sum((times - mean)**2)/real(nbuild))

      means(k) = mean
      stds(k) = std

      deallocate(times)
      deallocate(points)
    end do
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Save results to a file
    open(unit=10, file='../../../../benchmark/cosmokdtree_16_build_times.txt', status='replace')
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
