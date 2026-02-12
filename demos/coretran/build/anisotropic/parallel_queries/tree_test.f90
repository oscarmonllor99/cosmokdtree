program main
    use cosmokdtree
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

    !parameters
    integer :: ncpu
    integer :: j, i
    integer, parameter :: prec = 4 
    integer, parameter :: intkind = 4

    !input data
    real(kind=prec), allocatable :: L(:)
    real(kind=prec), allocatable :: x(:), y(:), z(:)
    real(kind=prec), allocatable :: points(:,:)
    integer(kind=intkind) :: n
    
    !tree results
    type(KDTreeNode), pointer :: root
    type(KDTreeResult) :: query
    integer(kind=intkind), allocatable :: indices(:)
    real(kind=prec), allocatable :: dist(:)

    !time
    integer*8 :: t1, t2, trate, tmax

    !queries
    real(kind=prec), allocatable :: ttarget(:,:)
    real(kind=prec) :: ball_radius
    integer :: k, nquery, ninside

    !file bigger than 2GB
    integer :: npart_save
    integer :: nsaves, nchecker
    character(len=1) :: ifile
    !!!!!!!!!!!!!!!!!!!!!!!!

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

    ! ! Read input data (FILE bigger than 2GB)
    ! n = 1000000000
    ! nsaves = 10
    ! nchecker = 0
    ! call system_clock(t1,trate,tmax)
    ! allocate(x(n), y(n), z(n))
    ! do i=1, nsaves
    !   write(ifile, '(I0)') i-1
    !   write(*,*) "Reading file:", trim('input/points_9_'//ifile//'.dat')
    !   open(unit=10, file=trim('input/points_9_'//ifile//'.dat'), form = 'unformatted')
    !   read(10) npart_save
    !   nchecker = nchecker + npart_save
    !   read(10) x((i-1)*npart_save+1:i*npart_save)
    !   read(10) y((i-1)*npart_save+1:i*npart_save)
    !   read(10) z((i-1)*npart_save+1:i*npart_save)
    !   close(10)
    ! enddo
    ! call system_clock(t2,trate,tmax)
    ! WRITE(*,*) "Time taken to read points:", float(t2 - t1)/float(trate), "seconds"
    ! WRITE(*,*) "Number of points:", n, nchecker
    ! WRITE(*,*) minval(x), maxval(x)
    ! WRITE(*,*) minval(y), maxval(y)
    ! WRITE(*,*) minval(z), maxval(z)

    ! Read input data
    call system_clock(t1,trate,tmax)
    open(unit=10, file='../input/points.dat', form = 'unformatted')
    read(10) n
    allocate(x(n), y(n), z(n))
    read(10) x
    read(10) y
    read(10) z
    close(10)
    call system_clock(t2,trate,tmax)
    ! WRITE(*,*) "Time taken to read points:", float(t2 - t1)/float(trate), "seconds"
    ! WRITE(*,*) "Number of points:", n

    allocate(points(n,3))
    points(:,1) = x
    points(:,2) = y
    points(:,3) = z

    ! Build the KD-Tree
    call system_clock(t1,trate,tmax)
    root => build_kdtree(points)
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to build KD-Tree:", float(t2 - t1)/float(trate), "seconds"

    ! Build the KD-Tree coretran
    x_coretran = x
    y_coretran = y
    z_coretran = z
    call system_clock(t1,trate,tmax)
    tree_coretran = KDTREE(x_coretran, y_coretran, z_coretran)
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "CORETRAN Time taken to build KD-Tree:", float(t2 - t1)/float(trate), "seconds"

    ! QUERIES
    allocate(L(3))
    L = 100.
    nquery = 100000
    allocate(ttarget(nquery,3))
    do j = 1, 3 
      do i = 1, nquery
        ttarget(i,j) = rand() * L(j)  ! Random target point
      end do
    end do

    ! knn-search 
    k = 100
    call system_clock(t1,trate,tmax)
    !$omp parallel shared(root, ttarget, k) private(i, query)
    !$omp do
    do i = 1, nquery
      query = knn_search(root, ttarget(i,:), k)
    end do
    !$omp end do
    !$omp end parallel
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to find nearests neighbors:", float(t2 - t1)/float(trate), "seconds"

    !CORETRAN KNN TEST
    call system_clock(t1,trate,tmax)
    !$omp parallel shared(tree_coretran, x_coretran, y_coretran, z_coretran, ttarget, k) private(i, DA)
    !$omp do
    do i = 1, nquery
      DA = SEARCH%KNEAREST(tree_coretran, x_coretran, y_coretran, z_coretran, &
              XQUERY=DBLE(ttarget(i,1)), YQUERY=DBLE(ttarget(i,2)),  &
              ZQUERY=DBLE(ttarget(i,3)), K=k)
    end do
    !$omp end do
    !$omp end parallel
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "CORETRAN Time taken to find nearests neighbors:", float(t2 - t1)/float(trate), "seconds"

    ! ball-search
    ball_radius = 1.
    call system_clock(t1,trate,tmax)
    !$omp parallel shared(root, ttarget, ball_radius) private(i, query)
    !$omp do
    do i = 1, nquery
      query = ball_search(root, ttarget, ball_radius)
    end do
    !$omp end do
    !$omp end parallel
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to find points in ball:", float(t2 - t1)/float(trate), "seconds"
    
    !CORETRAN BALL TEST
    call system_clock(t1,trate,tmax)
    !$omp parallel shared(tree_coretran, x_coretran, y_coretran, z_coretran, ttarget, ball_radius) private(i, DA)
    !$omp do
    do i=1, nquery
      DA = SEARCH%KNEAREST(tree_coretran, x_coretran, y_coretran, z_coretran, &
                XQUERY=DBLE(ttarget(i,1)), YQUERY=DBLE(ttarget(i,2)),  &
                ZQUERY=DBLE(ttarget(i,3)), RADIUS=DBLE(ball_radius))  
    end do
    !$omp end do
    !$omp end parallel
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "CORETRAN Time taken to find points within the ball:", float(t2 - t1)/float(trate), "seconds"

    ! Deallocate memory
    deallocate(x, y, z)
    deallocate(points)
    deallocate(root)
end program main
