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
    real(kind=prec), allocatable :: x(:), y(:), z(:)
    real(kind=prec), allocatable :: points(:,:)
    integer(kind=intkind) :: n, i, k
    integer(kind=intkind), allocatable :: ninputs(:)
    integer, allocatable :: exps(:)
    character(10) :: exp_str
    !file bigger than 2GB
    integer :: npart_save
    integer :: nsaves, nchecker
    character(len=1) :: ifile
    !!!!!!!!!!!!!!!!!!!!!!!!
    
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
    nbuild = 5
    allocate(L(ndim))
    allocate(ninputs(4))
    allocate(means(4))
    allocate(stds(4))
    allocate(exps(4))
    exps = [6, 7, 8, 9]  ! Corresponding to 10^6, 10^7, 10^8, 10^9
    ninputs(1) = 1000000
    ninputs(2) = 10000000
    ninputs(3) = 100000000
    ninputs(4) = 1000000000
    means = 0.0
    stds = 0.0

    !file bigger than 2GB
    nsaves = 10
    !!!!!!!!!!!!!!!!!!!!!!!
    
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    do k = 1, size(ninputs)
      n = ninputs(k)
      allocate(x(n), y(n), z(n))
      print *, "Building KD-Tree with ", n, " points."
      write(exp_str, '(I0)') exps(k)
      if (n .lt. 1000000000) then
        open(unit=10, file='../../../../benchmark/examples/points_'&
            // trim(adjustl(exp_str)) // '.dat', &
              form = 'unformatted')

        ! random uniform points in a cube of size L
        read(10) !n
        read(10) x
        read(10) y
        read(10) z
        close(10)
      else 
        nchecker = 0
        do i = 1, nsaves
          write(ifile, '(I0)') i-1
          open(unit=10, file=trim('../../../../benchmark/examples&
                             &/points_'// trim(adjustl(exp_str))//'_'//ifile//'.dat'), form = 'unformatted')
          read(10) npart_save
          nchecker = nchecker + npart_save
          read(10) x((i-1)*npart_save+1:i*npart_save)
          read(10) y((i-1)*npart_save+1:i*npart_save)
          read(10) z((i-1)*npart_save+1:i*npart_save)
          close(10)
        end do
      end if
    
      allocate(points(n,ndim))
      points(:,1) = x
      points(:,2) = y
      points(:,3) = z
      deallocate(x,y,z)

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
    open(unit=10, file='../../../../benchmark/inhomo_cosmokdtree_build_times.txt', status='replace')
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
