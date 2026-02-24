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

    !input data
    integer :: n, i, k, ndim
    real, allocatable :: x(:), y(:), z(:)
    integer, allocatable :: ninputs(:)
    integer, allocatable :: exps(:)
    character(10) :: exp_str
  
    !time
    integer*8 :: t1, t2, trate, tmax
    real*4, allocatable :: times(:), means(:), stds(:)
    real*4 :: mean, std
    integer ::  nbuild

    !file bigger than 2GB
    integer :: npart_save
    integer :: nsaves, nchecker
    character(len=1) :: ifile
    !!!!!!!!!!!!!!!!!!!!!!!!

    !coretran
    TYPE(KDTREE) :: tree_coretran
    REAL(R64), ALLOCATABLE :: xtree(:), ytree(:), ztree(:)

    ! Set the number of dimensions
    ndim = 3
    nbuild = 5
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

      call allocate(xtree, n)
      call allocate(ytree, n)
      call allocate(ztree, n)
      xtree = x
      ytree = y
      ztree = z

      ! Build the KD-Tree and time it
      allocate(times(nbuild))  
      do i = 1, nbuild
        write(*,*) " Build iteration ", i, " of ", nbuild
        call system_clock(t1,trate,tmax)
        tree_coretran = KDTREE(xtree, ytree, ztree)
        call system_clock(t2,trate,tmax)
        times(i) = float(t2 - t1)/float(trate)
        call tree_coretran%deallocate()
      end do
      
      mean = sum(times)/real(nbuild)
      std = sqrt(sum((times - mean)**2)/real(nbuild))

      means(k) = mean
      stds(k) = std
      
      deallocate(times)
      deallocate(x,y,z)
      deallocate(xtree, ytree, ztree)
    end do
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Save results to a file
    open(unit=10, file='../../../../benchmark/inhomo_coretran_build_times.txt', status='replace')
    write(10,*) 'N_points  Mean_time(s)  Std_dev(s)'
    do k = 1, size(ninputs)
      write(10,*) ninputs(k), means(k), stds(k)
    end do
    close(10)
    deallocate(ninputs)
    deallocate(means)
    deallocate(stds)  

end program main
