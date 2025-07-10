! Main friends of friends subroutine
module FoF_test
    implicit none
    contains
    subroutine friends_of_friends(x,y,z,n,ll,maxGroups,root,partGroup,Ngroup,NfriendGroup, &
                                  ngrid,dx,head,next)
        use cosmokdtree
        implicit none

        !input
        integer :: n, ngrid
        integer :: maxGroups
        real :: ll, dx
        real :: x(n), y(n), z(n)
        type(KDTreeNode), pointer :: root
        integer :: head(ngrid**3), next(n)

        !local
        integer :: ip, ip2, ip3
        logical :: already_looked(n), already_friended(n)
        integer :: friends(n)
        integer :: queue(n)
        integer :: Nfriend
        integer :: poin

        !results
        integer :: partGroup(n)
        integer :: NfriendGroup(maxGroups)
        integer :: Ngroup

        partGroup = 0
        Ngroup = 0
        NfriendGroup = 0

        !FoF local variables
        already_looked = .false.
        already_friended = .false.
        friends = 0
        Nfriend = 0
        do ip = 1, n

            if (already_looked(ip)) cycle
            already_looked(ip) = .true.

            ! CHOOSE SEARCH METHOD !!!!!!!!!
            !call search_kdtree(x(ip), y(ip), z(ip), n, ll, root, friends, Nfriend)
            !call search_grid(x,y,z,ip,n,ll,ngrid,dx,head,next,friends,Nfriend)
            call search_brute_force(x, y, z, ip, n, ll, friends, Nfriend)

            if (Nfriend .le. 1) then
                cycle
            end if

            ! new group is created
            Ngroup = Ngroup + 1
            if (Ngroup > maxGroups) then
                write(*,*) "STOP: Maximum number of groups exceeded", maxGroups
                stop
            end if

            do ip2 = 1, Nfriend
                NfriendGroup(Ngroup) = NfriendGroup(Ngroup) + 1
                partGroup(friends(ip2)) = Ngroup
                already_friended(friends(ip2)) = .true.
            end do

            !friends-of-friends loop: queue of friends to look for friends of friends
            poin = 0
            queue(1:Nfriend) = friends(1:Nfriend)
            do while (poin < NfriendGroup(Ngroup))
                poin = poin + 1
                ip2 = queue(poin)

                if (already_looked(ip2)) cycle
                
                ! CHOOSE SEARCH METHOD !!!!!!!!!
                !call search_kdtree(x(ip2), y(ip2), z(ip2), n, ll, root, friends, Nfriend)
                !call search_grid(x,y,z,ip2,n,ll,ngrid,dx,head,next,friends,Nfriend)
                call search_brute_force(x, y, z, ip2, n, ll, friends, Nfriend)

                do ip3 = 1, Nfriend
                    if (already_looked(friends(ip3))) cycle
                    if (already_friended(friends(ip3))) cycle
                    NfriendGroup(Ngroup) = NfriendGroup(Ngroup) + 1
                    partGroup(friends(ip3)) = Ngroup
                    already_friended(friends(ip3)) = .true.
                    !add to queue to look for friends of friends
                    queue(NfriendGroup(Ngroup)) = friends(ip3)
                end do

                already_looked(ip2) = .true.

            end do
        enddo
    end subroutine

    subroutine search_kdtree(x, y, z, n, R, root, friends, counter)
        use cosmokdtree
        implicit none
        !in
        real :: x, y, z, R
        type(KDTreeNode), pointer :: root
        integer :: n
        !local
        type(KDTreeResult) :: query
        real :: ttarget(3)
        !out
        integer :: friends(n)
        integer :: counter

        !query point
        ttarget(1) = x
        ttarget(2) = y
        ttarget(3) = z

        !query
        query = ball_search(root, ttarget, R)

        !get results
        counter = size(query%idx)
        friends(1:counter) = query%idx

    end subroutine
    
    subroutine search_grid(x, y, z, ip, n, R, ngrid, dx, head, next, friends, counter)
        implicit none
        !in
        real :: x(n), y(n), z(n)
        integer :: n, ip, ngrid
        real :: R, dx
        integer :: head(ngrid**3), next(n)
        !out
        integer :: friends(n)
        integer :: counter
        !local
        integer :: i, i1d, ix, jy, kz, ix2, jy2, kz2

        counter = 0
        friends = 0

        ! This particle's cell
        ix = int(x(ip) / dx) + 1
        jy = int(y(ip) / dx) + 1
        kz = int(z(ip) / dx) + 1

        do ix2 = ix-1,ix+1
            if (ix2 < 1 .or. ix2 > ngrid) cycle
            do jy2 = jy-1,jy+1
                if (jy2 < 1 .or. jy2 > ngrid) cycle
                do kz2 = kz-1,kz+1
                    if (kz2 < 1 .or. kz2 > ngrid) cycle

                    ! 3D index to 1D in fortran order
                    i1d = ix2 + (jy2 - 1)*ngrid + (kz2 - 1)*ngrid*ngrid

                    ! Recover all particles in this cell -> linked list
                    i = head(i1d)
                    do while (i .ne. -1)
                        if (sqrt((x(i) - x(ip))**2 + (y(i) - y(ip))**2 + (z(i) - z(ip))**2) .le. R) then
                            counter = counter + 1
                            friends(counter) = i
                        end if
                        i = next(i)
                    end do
                end do
            end do
        end do

    end subroutine
    
    subroutine search_brute_force(x, y, z, ip, n, R, friends, counter)
        implicit none
        !in
        real :: x(n), y(n), z(n)
        integer :: n, ip
        real :: R
        !out
        integer :: friends(n)
        integer :: counter
        !local
        integer :: i

        counter = 0
        friends = 0

        do i = 1, n 
            if (sqrt((x(i) - x(ip))**2 + (y(i) - y(ip))**2 + (z(i) - z(ip))**2) .le. R) then
                counter = counter + 1
                friends(counter) = i
            end if
        end do
        
    end subroutine
    
end module FoF_test



