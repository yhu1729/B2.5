program test_sdc_quad
!   Standalone check of b2mod_sdc_quad. Build and run with
!     gfortran -o test_sdc_quad ../../modules/b2mod_sdc_quad.F90 test_sdc_quad.F90 && ./test_sdc_quad
!   Exit status 0 iff all node counts pass.
    use b2mod_sdc_quad
    implicit none
    integer :: m
    logical :: ok, all_ok
    real(sdc_dp) :: q(SDC_MAX_NODES,SDC_MAX_NODES), qd(SDC_MAX_NODES,SDC_MAX_NODES)
    real(sdc_dp) :: c(SDC_MAX_NODES), qd_ref(2,2)
    integer :: ierr, i

    all_ok = .true.
    do m = 1, SDC_MAX_NODES
        call sdc_quad_selftest(m, 1.0e-13_sdc_dp, ok, report=.true.)
        all_ok = all_ok .and. ok
    end do
!   LU trick sanity: Q_delta lower triangular with positive diagonal
    do m = 1, SDC_MAX_NODES
        call sdc_radau_nodes(m, c, ierr)
        call sdc_collocation_matrix(m, c, q)
        call sdc_qdelta_lu(m, q, qd)
        do i = 1, m
            if (qd(i,i) <= 0.0_sdc_dp) all_ok = .false.
            if (i < m) then
                if (maxval(abs(qd(i,i+1:m))) > 1.0e-14_sdc_dp) all_ok = .false.
            end if
        end do
        print '(a,i0,a,4es22.14)', 'Q_delta diag M=', m, ': ', (qd(i,i), i=1,m)
    end do
!   closed forms: M = 1: Q = [1], Q_delta = [1].
!   M = 2 (Radau IIA, c = 1/3, 1): Q = [[5/12, -1/12], [3/4, 1/4]],
!   Q^T = L U with l21 = -1/5, U = [[5/12, 3/4], [0, 2/5]],
!   Q_delta = U^T = [[5/12, 0], [3/4, 2/5]].
    call sdc_radau_nodes(1, c, ierr)
    call sdc_collocation_matrix(1, c, q)
    call sdc_qdelta_lu(1, q, qd)
    if (abs(qd(1,1) - 1.0_sdc_dp) > 1.0e-14_sdc_dp) then
        all_ok = .false.
        print '(a,es12.3)', 'Q_delta M=1 mismatch', abs(qd(1,1) - 1.0_sdc_dp)
    end if
    call sdc_radau_nodes(2, c, ierr)
    call sdc_collocation_matrix(2, c, q)
    call sdc_qdelta_lu(2, q, qd)
    qd_ref(1,1) = 5.0_sdc_dp/12.0_sdc_dp
    qd_ref(1,2) = 0.0_sdc_dp
    qd_ref(2,1) = 0.75_sdc_dp
    qd_ref(2,2) = 0.4_sdc_dp
    if (maxval(abs(qd(1:2,1:2) - qd_ref)) > 1.0e-14_sdc_dp) then
        all_ok = .false.
        print '(a,es12.3)', 'Q_delta M=2 mismatch', maxval(abs(qd(1:2,1:2) - qd_ref))
    end if
    if (all_ok) then
        print '(a)', 'test_sdc_quad: PASS'
    else
        print '(a)', 'test_sdc_quad: FAIL'
        stop 1
    end if
end program test_sdc_quad
