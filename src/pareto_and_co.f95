module pareto_and_co
    use iso_c_binding
    
    implicit none
    
    private
    public pareto, dominate, dominated
    
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
        real(kind=c_double), dimension(nind, nobj) :: Xrow
        integer(kind=c_int), intent(out), dimension(nind) :: Ft
        integer(kind=c_int), dimension(nind) :: ipop
        logical, dimension(nind) :: is_efficient
        integer, allocatable :: ileft(:), is_dominated(:)
        integer :: i, j, k, nleft
        
        ipop = [(i, i=1, nind)]
        is_efficient = .TRUE.
        
        do i = 1, nind
            if(is_efficient(i)) then
              nleft = count(is_efficient)
              allocate(ileft(nleft))
              allocate(is_dominated(nleft))
              ileft = pack(ipop, is_efficient)
              call dominated(X(i,:), X(ileft,:), nleft, nobj, is_dominated)
              is_efficient(ileft) = is_dominated == 0 .AND. is_efficient(ileft)
              deallocate(ileft,is_dominated)
            end if           
        end do
        Ft = merge(1, 0, is_efficient)
    
    end subroutine pareto
    
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
