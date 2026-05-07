module pareto_and_co
    use iso_c_binding
    
    implicit none
    
    private
    public pareto, pareto_2d, pareto_3d, dominate, dominated
    
    contains
    
    subroutine dominate(matobj, nind, nobj, f) bind(C, name="dominate_")
        integer(kind = c_int), intent(in), value                :: nind, nobj
        real(kind = c_double), intent(in), dimension(nind,nobj) :: matobj
        integer(kind = c_int), intent(out), dimension(nind)     :: f
        integer, dimension(nind) :: ionion
        integer, allocatable :: ftmp(:), ipeel(:), ileft(:)
        logical, dimension(nind) :: feq0
        integer :: k, i, n, m
        
        ionion = [(i, i=1, nind)] ! overall onion index
        f = 0                     ! if 0 => need to be peeled
        k = 1                     ! peel counter
        feq0 = f == 0
        do while (any(feq0))
            n = count(feq0)
            allocate(ftmp(n))
            allocate(ileft(n))
            ileft = pack(ionion, feq0) ! the left onion
            call pareto(matobj(ileft,:), n, nobj, ftmp)
            m = sum(ftmp)
            allocate(ipeel(m))
            ipeel = pack(ionion(ileft), ftmp == 1) ! the peel 
            f(ipeel) = k 
            k = k + 1
            feq0 = f == 0
            deallocate(ftmp, ipeel, ileft)
        end do

    end subroutine dominate

    subroutine pareto(X, nind, nobj, Ft) bind(C, name="pareto_")
        integer(kind=c_int), intent(in), value :: nind, nobj
        real(kind=c_double), intent(in), dimension(nind, nobj) :: X
        integer(kind=c_int), intent(out), dimension(nind) :: Ft
    
        integer :: i, k, npf, kobj
        integer, allocatable :: pf(:)        ! indices of the current Pareto points
        real(kind=c_double) :: xi_k, xj_k
        logical :: dominated_i, strictly_better_found, equal_so_far
        logical :: any_strict_in_pf
    
        allocate(pf(nind))
        npf = 0
        Ft = 0
    
        do i = 1, nind
            dominated_i = .false.
            ! Comparison of X(i,:) with each Pareto point already retained pf(1..npf)
            do k = 1, npf
                ! Test: does pf(k) dominate i?
                !   pf(k) dominates i <=> X(pf(k),:) >= X(i,:) on all objectives
                !                        AND X(pf(k),:) > X(i,:) on at least one.
                ! As we maximize (cf. caRamel), "better" = larger.
                strictly_better_found = .false.
                equal_so_far = .true.
                any_strict_in_pf = .true.
                do kobj = 1, nobj
                    xi_k = X(i, kobj)
                    xj_k = X(pf(k), kobj)
                    if (xj_k < xi_k) then
                        ! pf(k) is worse on this objective => does not dominate i
                        any_strict_in_pf = .false.
                        exit
                    else if (xj_k > xi_k) then
                        strictly_better_found = .true.
                        equal_so_far = .false.
                    end if
                end do
                if (any_strict_in_pf) then
                    ! pf(k) >= i everywhere; dominates if strictly better somewhere
                    ! or if equal everywhere (duplicate: we keep the first one encountered,
                    ! so i is marked as dominated).
                    if (strictly_better_found .or. equal_so_far) then
                        dominated_i = .true.
                        exit
                    end if
                end if
            end do
    
            if (.not. dominated_i) then
                ! i is non-dominated so far. We add it to pf,
                ! but first we remove from pf the points that i now dominates.
                call prune_and_add(X, nind, nobj, i, pf, npf)
            end if
        end do
    
        Ft = 0
        do k = 1, npf
            Ft(pf(k)) = 1
        end do
        deallocate(pf)
    end subroutine pareto
    
    subroutine prune_and_add(X, nind, nobj, i, pf, npf)
        use iso_c_binding
        integer(kind=c_int), intent(in) :: nind, nobj, i
        real(kind=c_double), intent(in) :: X(nind, nobj)
        integer, intent(inout) :: pf(:)
        integer, intent(inout) :: npf
        integer :: k, w, kobj
        real(kind=c_double) :: xi_k, xj_k
        logical :: i_dominates_k, strictly_better, equal_so_far
    
        w = 0
        do k = 1, npf
            ! does i dominate pf(k)?
            i_dominates_k = .true.
            strictly_better = .false.
            equal_so_far = .true.
            do kobj = 1, nobj
                xi_k = X(i, kobj)
                xj_k = X(pf(k), kobj)
                if (xi_k < xj_k) then
                    i_dominates_k = .false.
                    exit
                else if (xi_k > xj_k) then
                    strictly_better = .true.
                    equal_so_far = .false.
                end if
            end do
            if (i_dominates_k .and. strictly_better) then
                ! pf(k) is dominated by i => we remove it (not rewritten into pf)
                cycle
            else
                w = w + 1
                pf(w) = pf(k)
            end if
        end do
        w = w + 1
        pf(w) = i
        npf = w
    end subroutine prune_and_add

    subroutine pareto_2d(X, nind, Ft) bind(C, name="pareto_2d_")
        integer(kind=c_int), intent(in), value :: nind
        real(kind=c_double), intent(in), dimension(nind, 2) :: X
        integer(kind=c_int), intent(out), dimension(nind) :: Ft
    
        integer :: i
        real(kind=c_double) :: f1_prev, f2_max
    
        if (nind <= 0) return
    
        Ft = 0
        ! The first point (f1 max) is necessarily non-dominated.
        Ft(1) = 1
        f1_prev = X(1, 1)
        f2_max  = X(1, 2)
    
        do i = 2, nind
            ! Thanks to the sorting (f1 decreasing, then f2 decreasing):
            !  - X(i,1) <= X(i-1,1)
            !  - if X(i,1) == f1_prev : f1 duplicate, X(i,2) <= f2_max 
            !    so dominated (or equal to an already retained point) => Ft(i)=0
            !  - if X(i,1) <  f1_prev :
            !      non-dominated if and only if X(i,2) > f2_max (strict)
            if (X(i, 1) < f1_prev) then
                if (X(i, 2) > f2_max) then
                    Ft(i) = 1
                    f2_max = X(i, 2)
                end if
                f1_prev = X(i, 1)
            end if
            ! otherwise: f1 equal to the previous one, 
            !      point dominated/equal => Ft(i) remains 0
        end do
    end subroutine pareto_2d

    subroutine pareto_3d(X, nind, Ft) bind(C, name="pareto_3d_")
        !-------------------------------------------------------------
        ! X is assumed lexicographically sorted DESCENDING on
        ! (f1, f2, f3). We maintain a staircase of current Pareto
        ! points, indexed by decreasing f2 (f3 strictly increasing).
        ! Storage: index array 'stair(:)' of size npf.
        !-------------------------------------------------------------
        integer(kind=c_int), intent(in), value :: nind
        real(kind=c_double), intent(in), dimension(nind, 3) :: X
        integer(kind=c_int), intent(out), dimension(nind) :: Ft
    
        integer, allocatable :: stair(:)
        integer :: i, npf, lo, hi, mid, k, j, new_npf
        real(kind=c_double) :: f2p, f3p
    
        Ft = 0
        if (nind <= 0) return
    
        allocate(stair(nind))
    
        ! First point: always non-dominated (largest f1)
        Ft(1) = 1
        stair(1) = 1
        npf = 1
    
        do i = 2, nind
            f2p = X(i, 2)
            f3p = X(i, 3)
    
            !----- Binary search: smallest k such that X(stair(k),2) < f2p -----
            ! stair is sorted by f2 DESCENDING, so we look for the boundary.
            lo = 1
            hi = npf + 1
            do while (lo < hi)
                mid = (lo + hi) / 2
                if (X(stair(mid), 2) < f2p) then
                    hi = mid
                else
                    lo = mid + 1
                end if
            end do
            k = lo ! 1..npf+1 ; indices < k have f2 >= f2p
    
            !----- Dominance test -----
            ! The potential dominator candidate is stair(k-1) (largest f3
            ! among those with f2 >= f2p, since f3 increases 
            ! with index in stair).
            if (k > 1) then
                if (X(stair(k-1), 3) >= f3p) then
                    ! Dominated (weakly) by stair(k-1).
                    ! Since stair(k-1) has f1 >= f1(i) (sorting) and f2 >= f2p
                    ! and f3 >= f3p, and (f1,f2,f3) cannot be strictly equal
                    ! (otherwise i would come after in descending lex order
                    ! and stair(k-1) already retained)
                    ! => i is strictly dominated.
                    cycle
                end if
            end if
    
            !----- i is non-dominated: insert and purge -----
            Ft(i) = 1
    
            ! Purge: all stair(k), stair(k+1), ... with f3 <= f3p
            ! are dominated by i (i has f2 > their f2, and f3 >= their f3).
            ! Since f3 increases with index, this is a contiguous prefix
            ! starting from k. We remove them by revoking Ft.
            j = k
            do while (j <= npf)
                if (X(stair(j), 3) <= f3p) then
                    Ft(stair(j)) = 0
                    j = j + 1
                else
                    exit
                end if
            end do
    
            ! Rebuild the stair array: elements [1..k-1], then i,
            ! then [j..npf] (those not purged).
            ! Done in place via shifting.
            if (j == k) then
                ! No purge: simple insertion at position k (right shift)
                ! Shift stair(k..npf) one slot to the right
                ! (explicit loop to stay in place)
                ! npf becomes npf+1
                new_npf = npf + 1
                ! Shift from end to start to avoid overwriting
                ! Loop from npf down to k:
                block
                    integer :: t
                    do t = npf, k, -1
                        stair(t+1) = stair(t)
                    end do
                end block
                stair(k) = i
                npf = new_npf
            else
                ! Purge (j - k) elements starting at k, then insert i at k.
                ! New size = npf - (j - k) + 1
                ! If j > k+1, shift stair(j..npf) to stair(k+1..)
                ! If j < k+1 impossible here. If j == k+1, it is exactly
                ! a replacement: stair(k) <- i.
                if (j == k + 1) then
                    stair(k) = i
                    ! npf unchanged
                else
                    ! Shift stair(j..npf) to stair(k+1..k+1+(npf-j))
                    block
                        integer :: t, shift
                        shift = j - (k + 1)   ! > 0
                        do t = j, npf
                            stair(t - shift) = stair(t)
                        end do
                    end block
                    stair(k) = i
                    npf = npf - (j - k) + 1
                end if
            end if
        end do
    
        deallocate(stair)
    end subroutine pareto_3d

    subroutine dominated(Xi, X, nind, nobj, is_dominated) bind(C, name="dominated_")
        
        integer(kind=c_int), intent(in), value :: nind, nobj
        real(kind=c_double), intent(in), dimension(nobj) :: Xi
        real(kind=c_double), intent(in), dimension(nind, nobj) :: X
        integer(kind=c_int), intent(out), dimension(nind) :: is_dominated
        logical, dimension(nind) :: is_efficient
        integer :: j
        
        do j = 1, nind
          is_efficient(j) = any(X(j,:) > Xi) .OR. all(X(j,:) == Xi)
        end do
        is_dominated = merge(0, 1, is_efficient)
        
    end subroutine dominated

end module
