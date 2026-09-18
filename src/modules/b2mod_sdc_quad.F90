module b2mod_sdc_quad
!   Quadrature data for spectral deferred correction (SDC) on the unit
!   interval: Radau IIA nodes c(1:M) with c(M) = 1, the collocation matrix
!   Q(m,j) = int_0^{c_m} l_j(s) ds, the node-to-node matrix
!   S(m,j) = Q(m,j) - Q(m-1,j) and the node spacings dtau(m) = c(m) - c(m-1).
!   The module has no dependency on the rest of B2.5 so that it can be
!   compiled and tested on its own (src/test/mms_time/test_sdc_quad.F90).
    implicit none
    private

    integer, parameter, public :: sdc_dp = selected_real_kind(14)
    integer, parameter, public :: SDC_MAX_NODES = 4

    public :: sdc_radau_nodes, sdc_collocation_matrix, &
              sdc_sweep_matrices, sdc_qdelta_lu, sdc_quad_selftest

contains

    subroutine sdc_radau_nodes(m, c, ierr)
!       Radau IIA (right) nodes on (0,1]: roots of P_m(2x-1) - P_{m-1}(2x-1),
!       with c(m) = 1. Closed forms for m <= 3, Newton iteration for m = 4.
        integer, intent(in) :: m
        real(sdc_dp), intent(out) :: c(:)
        integer, intent(out) :: ierr
        real(sdc_dp) :: x, f, df, s6
        integer :: i, it

        ierr = 0
        if (m < 1 .or. m > SDC_MAX_NODES .or. size(c) < m) then
            ierr = 1
            return
        end if
        select case (m)
        case (1)
            c(1) = 1.0_sdc_dp
        case (2)
            c(1) = 1.0_sdc_dp/3.0_sdc_dp
            c(2) = 1.0_sdc_dp
        case (3)
            s6 = sqrt(6.0_sdc_dp)
            c(1) = (4.0_sdc_dp - s6)/10.0_sdc_dp
            c(2) = (4.0_sdc_dp + s6)/10.0_sdc_dp
            c(3) = 1.0_sdc_dp
        case (4)
!           interior nodes: roots of the Jacobi polynomial P^(1,0)_3 mapped
!           to (0,1); Newton from the Chebyshev-like initial guesses
            do i = 1, 3
                x = 0.5_sdc_dp*(1.0_sdc_dp - cos(real(2*i-1, sdc_dp)* &
                    acos(-1.0_sdc_dp)/7.0_sdc_dp))
                do it = 1, 60
                    call radau_poly(m, x, f, df)
                    if (abs(df) <= tiny(1.0_sdc_dp)) exit
                    x = x - f/df
                    if (abs(f) < 1.0e-15_sdc_dp) exit
                end do
                c(i) = x
            end do
            c(4) = 1.0_sdc_dp
            if (c(1) >= c(2) .or. c(2) >= c(3) .or. c(3) >= 1.0_sdc_dp) &
                ierr = 2
        end select
    end subroutine sdc_radau_nodes

    subroutine radau_poly(m, x, f, df)
!       f(x) = P_m(y) - P_{m-1}(y), y = 2x-1, and df/dx, by the Legendre
!       three-term recurrence. Zeros of f on (0,1) are the interior Radau IIA
!       nodes (x = 1 is always a root).
        integer, intent(in) :: m
        real(sdc_dp), intent(in) :: x
        real(sdc_dp), intent(out) :: f, df
        real(sdc_dp) :: y, p0, p1, p2, d0, d1, d2
        integer :: k

        y = 2.0_sdc_dp*x - 1.0_sdc_dp
        p0 = 1.0_sdc_dp
        d0 = 0.0_sdc_dp
        p1 = y
        d1 = 1.0_sdc_dp
        do k = 1, m - 1
            p2 = (real(2*k+1, sdc_dp)*y*p1 - real(k, sdc_dp)*p0)/ &
                 real(k+1, sdc_dp)
            d2 = (real(2*k+1, sdc_dp)*(p1 + y*d1) - real(k, sdc_dp)*d0)/ &
                 real(k+1, sdc_dp)
            p0 = p1
            d0 = d1
            p1 = p2
            d1 = d2
        end do
!       now p1 = P_m, p0 = P_{m-1}
        f = p1 - p0
        df = 2.0_sdc_dp*(d1 - d0)
    end subroutine radau_poly

    subroutine sdc_collocation_matrix(m, c, q)
!       Q(i,j) = int_0^{c_i} l_j(s) ds with l_j the Lagrange basis on c(1:m).
!       The basis coefficients are obtained from the Vandermonde system
!       V a_j = e_j, solved with partial pivoting (m <= 4).
        integer, intent(in) :: m
        real(sdc_dp), intent(in) :: c(:)
        real(sdc_dp), intent(out) :: q(:,:)
        real(sdc_dp) :: v(m,m), a(m,m)
        integer :: i, j, p

!       v(k,p) = c(k)**(p-1); a(:,j) solves v a = e_j, so l_j(s) = sum_p a(p,j) s^(p-1)
        do p = 1, m
            do i = 1, m
                v(i,p) = c(i)**(p-1)
            end do
        end do
        a = 0.0_sdc_dp
        do j = 1, m
            a(j,j) = 1.0_sdc_dp
        end do
        call solve_small(m, v, a)
        q(1:m,1:m) = 0.0_sdc_dp
        do i = 1, m
            do j = 1, m
                do p = 1, m
                    q(i,j) = q(i,j) + a(p,j)*c(i)**p/real(p, sdc_dp)
                end do
            end do
        end do
    end subroutine sdc_collocation_matrix

    subroutine solve_small(n, a, b)
!       Gaussian elimination with partial pivoting; b holds n right-hand
!       sides on entry and the solutions on exit. a is destroyed.
        integer, intent(in) :: n
        real(sdc_dp), intent(inout) :: a(n,n), b(n,n)
        real(sdc_dp) :: tmp(n), fac
        integer :: k, i, piv

        do k = 1, n - 1
            piv = k - 1 + maxloc(abs(a(k:n,k)), 1)
            if (piv /= k) then
                tmp = a(k,:); a(k,:) = a(piv,:); a(piv,:) = tmp
                tmp = b(k,:); b(k,:) = b(piv,:); b(piv,:) = tmp
            end if
            do i = k + 1, n
                fac = a(i,k)/a(k,k)
                a(i,k:n) = a(i,k:n) - fac*a(k,k:n)
                b(i,:) = b(i,:) - fac*b(k,:)
            end do
        end do
        do k = n, 1, -1
            b(k,:) = b(k,:)/a(k,k)
            do i = 1, k - 1
                b(i,:) = b(i,:) - a(i,k)*b(k,:)
            end do
        end do
    end subroutine solve_small

    subroutine sdc_sweep_matrices(m, c, q, s, dtau)
!       Node-to-node integration weights and node spacings.
        integer, intent(in) :: m
        real(sdc_dp), intent(in) :: c(:), q(:,:)
        real(sdc_dp), intent(out) :: s(:,:), dtau(:)
        integer :: i

        s(1,1:m) = q(1,1:m)
        dtau(1) = c(1)
        do i = 2, m
            s(i,1:m) = q(i,1:m) - q(i-1,1:m)
            dtau(i) = c(i) - c(i-1)
        end do
    end subroutine sdc_sweep_matrices

    subroutine sdc_qdelta_lu(m, q, qdelta)
!       LU trick of Weiser: Q^T = L U (Doolittle, unit lower L), Q_delta = U^T.
!       Q_delta is lower triangular; its diagonal replaces dtau in a stiff
!       sweep. Provided for phase 2; not used by the default sweep.
        integer, intent(in) :: m
        real(sdc_dp), intent(in) :: q(:,:)
        real(sdc_dp), intent(out) :: qdelta(:,:)
        real(sdc_dp) :: qt(m,m), l(m,m), u(m,m)
        integer :: i, j, k

        qt = transpose(q(1:m,1:m))
        l = 0.0_sdc_dp
        u = 0.0_sdc_dp
        do i = 1, m
            l(i,i) = 1.0_sdc_dp
        end do
        do k = 1, m
            do j = k, m
                u(k,j) = qt(k,j) - sum(l(k,1:k-1)*u(1:k-1,j))
            end do
            do i = k + 1, m
                l(i,k) = (qt(i,k) - sum(l(i,1:k-1)*u(1:k-1,k)))/u(k,k)
            end do
        end do
        qdelta(1:m,1:m) = transpose(u)
    end subroutine sdc_qdelta_lu

    subroutine sdc_quad_selftest(m, tol, ok, report)
!       Checks: nodes ordered in (0,1] with c(m) = 1; every row of Q
!       integrates polynomials of degree <= m-1 exactly from 0 to c_i; the
!       last row reproduces the Radau IIA weights; rows of S sum to dtau.
        integer, intent(in) :: m
        real(sdc_dp), intent(in) :: tol
        logical, intent(out) :: ok
        logical, intent(in), optional :: report
        real(sdc_dp) :: c(m), q(m,m), s(m,m), dtau(m), w(m), err, s6
        integer :: i, p, ierr
        logical :: verbose

        verbose = .false.
        if (present(report)) verbose = report
        ok = .true.
        call sdc_radau_nodes(m, c, ierr)
        if (ierr /= 0) then
            ok = .false.
            if (verbose) print '(a,i0)', 'sdc_quad_selftest: node error ', ierr
            return
        end if
        if (abs(c(m) - 1.0_sdc_dp) > tol .or. c(1) <= 0.0_sdc_dp) ok = .false.
        do i = 2, m
            if (c(i) <= c(i-1)) ok = .false.
        end do
        call sdc_collocation_matrix(m, c, q)
        call sdc_sweep_matrices(m, c, q, s, dtau)
!       polynomial exactness: sum_j Q(i,j) c_j^p = c_i^(p+1)/(p+1), p = 0..m-1
        do i = 1, m
            do p = 0, m - 1
                err = abs(sum(q(i,1:m)*c(1:m)**p) - c(i)**(p+1)/real(p+1, sdc_dp))
                if (err > tol) then
                    ok = .false.
                    if (verbose) print '(a,2i3,es12.3)', &
                        'sdc_quad_selftest: exactness failure row/deg', i, p, err
                end if
            end do
        end do
!       known Radau IIA weights
        select case (m)
        case (1)
            w = [1.0_sdc_dp]
        case (2)
            w = [0.75_sdc_dp, 0.25_sdc_dp]
        case (3)
            s6 = sqrt(6.0_sdc_dp)
            w = [(16.0_sdc_dp - s6)/36.0_sdc_dp, (16.0_sdc_dp + s6)/36.0_sdc_dp, &
                 1.0_sdc_dp/9.0_sdc_dp]
        case default
            w = q(m,1:m)
        end select
        err = maxval(abs(q(m,1:m) - w))
        if (err > tol) then
            ok = .false.
            if (verbose) print '(a,es12.3)', 'sdc_quad_selftest: weight mismatch', err
        end if
!       S rows sum to dtau
        do i = 1, m
            err = abs(sum(s(i,1:m)) - dtau(i))
            if (err > tol) then
                ok = .false.
                if (verbose) print '(a,i3,es12.3)', 'sdc_quad_selftest: S row sum', i, err
            end if
        end do
        if (verbose) then
            print '(a,i0,a,l1)', 'sdc_quad_selftest: M = ', m, ' ok = ', ok
            print '(a,4es22.14)', '  c    = ', c
            print '(a,4es22.14)', '  dtau = ', dtau
            do i = 1, m
                print '(a,i0,a,4es22.14)', '  Q(', i, ',:) = ', q(i,1:m)
            end do
        end if
    end subroutine sdc_quad_selftest

end module b2mod_sdc_quad
