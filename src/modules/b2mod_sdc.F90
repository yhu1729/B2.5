module b2mod_sdc
!   Spectral deferred correction (SDC) time integration wrapped around the
!   implicit B2.5 step b2mndt, and the driver-facing entry point b2step.
!
!   b2step has the argument list of b2mndt. With b2mndt_sdc = 0 it calls
!   b2mndt once (adding only the analytic MMS bookkeeping when
!   b2mndr_use_mms = 3). With b2mndt_sdc = 1 one time step [t_n, t_n + dt]
!   is integrated with Radau IIA nodes c(1:M), c(M) = 1, and K sweeps:
!
!     sweep 0: backward Euler node to node, dtim := dtau_m dt, psnl := U_{m-1}
!     sweep k: V (U_m - U_{m-1})/(dtau_m dt) = F(U_m) + b_m,
!              b_m = (1/dtau_m) sum_j S_mj F^k_j - F^k_m
!
!   With b2mndt_sdc_qdelta = 1 the sweep is preconditioned with the LU-trick
!   matrix Q_delta of Weiser (lower triangular, qd) instead of the node
!   spacings. Every node then steps from U_n with dtim := qd(m,m) dt:
!
!     sweep k: V (U_m - U_n)/(qd(m,m) dt) = F(U_m) + b_m,
!              b_m = (1/qd(m,m)) [ sum_j Q_mj F^k_j
!                                  + sum_{j<m} qd(m,j) (F^{k+1}_j - F^k_j) ]
!                    - F^k_m
!
!   with F^0 := 0 in sweep 0 (Q_delta predictor). This is exact for the
!   conserved time terms of b2npco/b2npmo/b2npht, which are differences
!   V (u - psnl)/dtim of na, m na ua and 3/2 n T. The qdelta = 0 path is
!   the node-to-node form above and is left unchanged.
!
!   F is the volume-integrated spatial residual without time term
!   (dv%resco0, resmo0, reshe0, reshi0) as stored by the Crank-Nicolson
!   branches of b2npco/b2npmo/b2npht when switch%b2mndt_ckn = 2, and b_m is
!   passed to those branches through psnl%res*0. The input b2mndt_ckn must
!   be 1 so that ts_factor stays 1; b2mndt_ckn is set to 2 only around the
!   internal b2mndt calls. F(U) after a node solve is obtained from an
!   evaluation-only b2mndt call (b2mndt_dummy = 2, rxf = 0, nstg = 1) whose
!   state changes are undone from a snapshot. The potential is algebraic:
!   it is solved in every node sweep and is not part of the quadrature.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use b2mod_types, only: R8
    use b2mod_switches, only: switches
    use b2us_geo, only: geometry
    use b2us_map, only: mapping
    use b2us_plasma, only: B2State, B2StateExt, B2Average, B2PlasmaSnapshot, &
                           createB2PlasmaSnapshot, getB2PlasmaSnapshot, &
                           putB2PlasmaSnapshot
    use b2mod_ad, only: ncall_b2mndt, ncall_b2news_, b2mndt_itcnt
    use b2mod_numerics_namelist, only: ts_factor, dtco, dtmo, dtee, dtei, &
                                       time_factor
    use b2mod_time, only: tim
    use b2mod_equation_sources, only: art_sna, art_smo, art_she, art_shi, &
                                      art_shn, art_sch
    use b2mod_boundary_sources, only: bv_na, bv_ua, bv_te, bv_ti, bv_tn, bv_po
    use b2mod_sdc_quad
    use b2mod_mms_analytic
    implicit none
    private
    public :: b2step

!   switches (b2mndt_sdc_*), read once
    integer, save :: sdc_on = 0
    integer, save :: sdc_nodes = 2
    integer, save :: sdc_sweeps = 2
    integer, save :: sdc_qdelta = 0
    integer, save :: sdc_nstg2 = 0
    integer, save :: sdc_iout = 0
    integer, save :: sdc_debug_stop = 0   ! b2mndt_sdc_debug_stop: abort after the first MMS source evaluation
    real(R8), save :: sdc_tol = 0.0_R8
    integer, save :: input_ckn = 1, input_bdf = 1

!   quadrature
    real(R8), save :: c(SDC_MAX_NODES), dtau(SDC_MAX_NODES)
    real(R8), save :: q(SDC_MAX_NODES,SDC_MAX_NODES)
    real(R8), save :: s(SDC_MAX_NODES,SDC_MAX_NODES)
    real(R8), save :: qd(SDC_MAX_NODES,SDC_MAX_NODES)   ! LU-trick Q_delta (sdc_qdelta = 1)

!   node storage: node(0) = U_n, node(m) = U_m with F(U_m) in res*0
    type(B2PlasmaSnapshot), allocatable, save :: node(:)
    type(B2PlasmaSnapshot), save :: work, psnl_entry, psnl_eval_save
!   F^k of the previous sweep (packed per equation), interior + guard cells
    real(R8), allocatable, save :: fk_co(:,:,:), fk_mo(:,:,:)
    real(R8), allocatable, save :: fk_he(:,:), fk_hi(:,:)
!   per-node MMS sources and boundary values (use_mms = 3)
    real(R8), allocatable, save :: nd_sna(:,:,:,:), nd_smo(:,:,:,:)
    real(R8), allocatable, save :: nd_she(:,:,:), nd_shi(:,:,:), nd_sch(:,:,:)
    real(R8), allocatable, save :: nd_bvna(:,:,:), nd_bvua(:,:,:)
    real(R8), allocatable, save :: nd_bvte(:,:), nd_bvti(:,:), nd_bvpo(:,:)

    logical, save :: initialised = .false.
    integer, save :: n_eval_total = 0
    integer, save :: n_step = 0
    real(R8), save :: last_sdc_res = 0.0_R8
    integer, save :: last_sdc_sweeps = 0

    external :: b2mndt, ipgeti, ipgetr, xertst, xerrab

contains

!-----------------------------------------------------------------------
    subroutine b2step(nout, nCv, nFc, nVx, ns, ismain, ismain0, nscx, &
                      nscxmax, iscx, itim, dtim, ntim, switch, geo, mpg, &
                      st, st_ext, st_avg, ierr)
        integer :: nout(0:*)
        integer, intent(in) :: nCv, nFc, nVx, ns, ismain, ismain0, nscx, &
                               nscxmax, itim, ntim
        integer :: iscx(0:*)
        real(R8), intent(in) :: dtim
        type(switches), intent(inout) :: switch
        type(geometry), intent(in) :: geo
        type(mapping), intent(inout) :: mpg
        type(B2State), intent(inout) :: st
        type(B2StateExt), intent(inout) :: st_ext
        type(B2Average), intent(inout) :: st_avg
        integer, intent(out) :: ierr

        integer :: m, k, ncall_entry, nCi
        logical :: need_eval, need_eval_node, converged
        real(R8) :: res

        ierr = 0
        nCi = mpg%nCi
        if (.not. initialised) call sdc_init(nCv, nFc, ns, switch, geo, mpg)

!       driver start-up call for Crank-Nicolson / BDF2-with-CN-start:
!       provide the sources at t_n, then let b2mndt do its dummy sweep
        if (switch%b2mndt_dummy == 1) then
            if (switch%use_mms == 3) &
                call mms_update_sources(tim, nout, nCv, nFc, nVx, ns, ismain, &
                    ismain0, nscx, nscxmax, iscx, itim, dtim, ntim, switch, &
                    geo, mpg, st, st_ext, st_avg)
            call b2mndt(nout, nCv, nFc, nVx, ns, ismain, ismain0, nscx, &
                        nscxmax, iscx, itim, dtim, ntim, switch, geo, mpg, &
                        st, st_ext, st_avg, ierr)
            return
        end if

        if (sdc_on == 0) then
!           plain implicit Euler / CN / BDF2 step
            if (switch%use_mms == 3) &
                call mms_update_sources(tim + dtim, nout, nCv, nFc, nVx, ns, &
                    ismain, ismain0, nscx, nscxmax, iscx, itim, dtim, ntim, &
                    switch, geo, mpg, st, st_ext, st_avg)
            call b2mndt(nout, nCv, nFc, nVx, ns, ismain, ismain0, nscx, &
                        nscxmax, iscx, itim, dtim, ntim, switch, geo, mpg, &
                        st, st_ext, st_avg, ierr)
            if (ierr /= 0) return
            n_step = n_step + 1
            if (switch%use_mms == 3) &
                call mms_write_error(tim + dtim, dtim, nCv, ns, geo, mpg, &
                    st%pl, b2mndt_itcnt, n_eval_total, 0.0_R8, 0)
            return
        end if

!       ---------------------------- SDC step ----------------------------
        ncall_entry = ncall_b2mndt
        node(0) = st%psnc
        psnl_entry = st%psnl

!       manufactured sources at the node times
        if (switch%use_mms == 3) then
            do m = 1, sdc_nodes
                call mms_update_sources(tim + c(m)*dtim, nout, nCv, nFc, nVx, &
                    ns, ismain, ismain0, nscx, nscxmax, iscx, itim, dtim, &
                    ntim, switch, geo, mpg, st, st_ext, st_avg)
                nd_sna(:,:,:,m) = art_sna
                nd_smo(:,:,:,m) = art_smo
                nd_she(:,:,m) = art_she
                nd_shi(:,:,m) = art_shi
                nd_sch(:,:,m) = art_sch
                nd_bvna(:,:,m) = bv_na
                nd_bvua(:,:,m) = bv_ua
                nd_bvte(:,m) = bv_te
                nd_bvti(:,m) = bv_ti
                nd_bvpo(:,m) = bv_po
            end do
        end if

        converged = .false.
        last_sdc_sweeps = 0
        last_sdc_res = 0.0_R8
        sweeps: do k = 0, sdc_sweeps
!           F^k from the previous sweep; not needed in sweep 0
            if (k > 0) then
                do m = 1, sdc_nodes
                    fk_co(:,:,m) = node(m)%resco0
                    fk_mo(:,:,m) = node(m)%resmo0
                    fk_he(:,m) = node(m)%reshe0
                    fk_hi(:,m) = node(m)%reshi0
                end do
            end if
            need_eval = k < sdc_sweeps .or. sdc_tol > 0.0_R8 .or. sdc_iout > 0
            do m = 1, sdc_nodes
                if (switch%use_mms == 3) call load_node_sources(m)
                call node_solve(m, k, nout, nCv, nFc, nVx, ns, ismain, ismain0, &
                    nscx, nscxmax, iscx, itim, dtim, ntim, switch, geo, mpg, &
                    st, st_ext, st_avg, ierr)
                if (ierr /= 0) then
!                   node solve produced an invalid state: roll back to U_n
!                   and let the driver retry the step with dtim/10
                    call putB2PlasmaSnapshot(st%pl, st%dv, node(0))
                    st%psnl = psnl_entry
                    ncall_b2mndt = ncall_entry + 1
                    write(*,'(a,i3,a,i3,a)') 'b2step: SDC node ', m, ' sweep ', k, &
                        ' failed, requesting timestep reduction'
                    ierr = 10
                    return
                end if
!               the Q_delta sweep needs F^{k+1}_j, j < m, within the sweep
                need_eval_node = need_eval .or. &
                                 (sdc_qdelta == 1 .and. m < sdc_nodes)
                if (need_eval_node) then
                    call eval_residual(nout, nCv, nFc, nVx, ns, ismain, ismain0, &
                        nscx, nscxmax, iscx, itim, dtim, ntim, switch, geo, mpg, &
                        st, st_ext, st_avg, node(m))
                end if
            end do
            last_sdc_sweeps = k + 1
            if (need_eval) then
                res = sdc_residual(ns, nCi, dtim, geo)
                last_sdc_res = res
                if (sdc_iout > 0) write(*,'(a,i9,a,i3,a,es12.4)') &
                    'b2step: SDC itim ', itim, ' sweep ', k, ' residual ', res
                if (sdc_tol > 0.0_R8 .and. res < sdc_tol) then
                    converged = .true.
                    exit sweeps
                end if
            end if
        end do sweeps
        if (sdc_iout > 0 .and. sdc_tol > 0.0_R8 .and. .not. converged) &
            write(*,'(a,i9,a,es12.4)') 'b2step: SDC itim ', itim, &
                ' not converged, residual ', last_sdc_res

!       pl already holds node M = U_{n+1}; restore the bookkeeping of the driver
        st%psnl = psnl_entry
        ncall_b2mndt = ncall_entry + 1
        n_step = n_step + 1
        ierr = 0
        if (switch%use_mms == 3) &
            call mms_write_error(tim + dtim, dtim, nCv, ns, geo, mpg, st%pl, &
                b2mndt_itcnt, n_eval_total, last_sdc_res, last_sdc_sweeps)
    end subroutine b2step

!-----------------------------------------------------------------------
    subroutine sdc_init(nCv, nFc, ns, switch, geo, mpg)
        integer, intent(in) :: nCv, nFc, ns
        type(switches), intent(inout) :: switch
        type(geometry), intent(in) :: geo
        type(mapping), intent(in) :: mpg
        integer :: m, ierr, av_batch_all
        real(R8) :: delta_min, delta_max
        logical :: ok

        input_ckn = switch%b2mndt_ckn
        input_bdf = switch%b2mndt_bdf
        call ipgeti('b2mndt_sdc', sdc_on)
        call ipgeti('b2mndt_sdc_nodes', sdc_nodes)
        call ipgeti('b2mndt_sdc_sweeps', sdc_sweeps)
        call ipgeti('b2mndt_sdc_qdelta', sdc_qdelta)
        call ipgeti('b2mndt_sdc_nstg2', sdc_nstg2)
        call ipgeti('b2mndt_sdc_iout', sdc_iout)
        call ipgeti('b2mndt_sdc_debug_stop', sdc_debug_stop)
        call ipgetr('b2mndt_sdc_tol', sdc_tol)
        call xertst(sdc_on == 0 .or. sdc_on == 1, 'b2mndt_sdc should be 0 or 1')
        call xertst(sdc_nodes >= 1 .and. sdc_nodes <= SDC_MAX_NODES, &
                    'b2mndt_sdc_nodes should be in 1..4')
        call xertst(sdc_sweeps >= 0, 'b2mndt_sdc_sweeps should be nonnegative')
        call xertst(sdc_qdelta == 0 .or. sdc_qdelta == 1, &
                    'b2mndt_sdc_qdelta should be 0 (node spacing) or 1 (LU trick)')
        call xertst(sdc_nstg2 >= 0, 'b2mndt_sdc_nstg2 should be nonnegative')
        call xertst(sdc_iout == 0 .or. sdc_iout == 1, 'b2mndt_sdc_iout should be 0 or 1')
        call xertst(sdc_tol >= 0.0_R8, 'b2mndt_sdc_tol should be nonnegative')

        if (switch%use_mms == 3) then
            call mms_init(nCv, ns, geo, mpg)
            call xertst(switch%equation_sources /= 0 .and. switch%boundary_sources /= 0, &
                'use_mms = 3 needs b2mndr_equation_sources = 1 and b2mndr_boundary_sources = 1')
            call xertst(switch%pot_eq == 1, 'use_mms = 3 needs b2news_poteq = 1')
            call xertst(switch%b2nppo_restr_po == 0.0_R8, &
                'use_mms = 3 needs b2nppo_restr_po = 0')
            call xertst(switch%b2mndt_style == 2, 'use_mms = 3 needs b2mndt_style = 2')
            call xertst(switch%tn_style == 0, 'use_mms = 3 needs b2mn_tn_style = 0')
        end if

        if (sdc_on == 1) then
            call xertst(switch%b2mndt_style == 2, 'b2mndt_sdc needs b2mndt_style = 2')
            call xertst(switch%b2mndt_ckn == 1 .and. switch%b2mndt_bdf == 1, &
                'b2mndt_sdc needs b2mndt_ckn = 1 and b2mndt_bdf = 1')
            call xertst(switch%b2mndt_dummy == 0, 'b2mndt_sdc needs b2mndt_dummy = 0')
            call xertst(switch%use_eirene == 0, 'b2mndt_sdc needs b2mndr_eirene = 0')
            call xertst(switch%density_control == 0, &
                'b2mndt_sdc needs b2mndt_density_control = 0')
            call xertst(switch%iav_run == 0, 'b2mndt_sdc needs b2mndt_av = 0')
            av_batch_all = 0
            call ipgeti('b2mndt_av_batch_all', av_batch_all)
            call xertst(av_batch_all == 0, 'b2mndt_sdc needs b2mndt_av_batch_all = 0')
            call xertst(switch%tn_style == 0, 'b2mndt_sdc needs b2mn_tn_style = 0')
            call xertst(switch%solve_keps == 0, 'b2mndt_sdc needs b2mndr_solve_keps = 0')
            call xertst(switch%BoRiS == 0.0_R8, 'b2mndt_sdc needs b2news_BoRiS = 0')
            call xertst(switch%facdrift_start == switch%facdrift_target .and. &
                        switch%facExB_start == switch%facExB_target .and. &
                        switch%facvis_start == switch%facvis_target, &
                'b2mndt_sdc needs stationary drift/ExB/viscosity ramps (start = target)')
            call xertst(all(dtco(0:ns-1,0:mpg%nnreg(0)) == 1.0_R8) .and. &
                        all(dtmo(0:ns-1,0:mpg%nnreg(0)) == 1.0_R8) .and. &
                        all(dtee(0:mpg%nnreg(0)) == 1.0_R8) .and. &
                        all(dtei(0:mpg%nnreg(0)) == 1.0_R8), &
                'b2mndt_sdc needs dtco = dtmo = dtee = dtei = 1')
            if (allocated(time_factor)) call xertst(all(time_factor == 1.0_R8), &
                'b2mndt_sdc needs time_factor = 1 (core_dt_suppression = 1)')
            delta_min = 0.0_R8
            delta_max = 0.0_R8
            call ipgetr('b2mndr_delta_min', delta_min)
            call ipgetr('b2mndr_delta_max', delta_max)
            call xertst(delta_min == 0.0_R8 .or. delta_max == 0.0_R8, &
                'b2mndt_sdc: disable the driver dt control (b2mndr_delta_min/max)')

            call sdc_radau_nodes(sdc_nodes, c, ierr)
            call xertst(ierr == 0, 'b2mndt_sdc: Radau node construction failed')
            call sdc_collocation_matrix(sdc_nodes, c, q)
            call sdc_sweep_matrices(sdc_nodes, c, q, s, dtau)
            call sdc_quad_selftest(sdc_nodes, 1.0e-12_R8, ok)
            call xertst(ok, 'b2mndt_sdc: quadrature self-test failed')
            qd = 0.0_R8
            if (sdc_qdelta == 1) then
                call sdc_qdelta_lu(sdc_nodes, q, qd)
                do m = 1, sdc_nodes
                    call xertst(qd(m,m) > 0.0_R8, &
                        'b2mndt_sdc: Q_delta has a nonpositive diagonal entry')
                end do
            end if

            allocate(node(0:sdc_nodes))
            do m = 0, sdc_nodes
                call createB2PlasmaSnapshot(nCv, nFc, ns, node(m))
            end do
            allocate(fk_co(nCv,0:ns-1,sdc_nodes), fk_mo(nCv,0:ns-1,sdc_nodes))
            allocate(fk_he(nCv,sdc_nodes), fk_hi(nCv,sdc_nodes))
            if (switch%use_mms == 3) then
                allocate(nd_sna(nCv,0:1,0:ns-1,sdc_nodes), nd_smo(nCv,0:3,0:ns-1,sdc_nodes))
                allocate(nd_she(nCv,0:3,sdc_nodes), nd_shi(nCv,0:3,sdc_nodes), &
                         nd_sch(nCv,0:3,sdc_nodes))
                allocate(nd_bvna(nCv,0:ns-1,sdc_nodes), nd_bvua(nCv,0:ns-1,sdc_nodes))
                allocate(nd_bvte(nCv,sdc_nodes), nd_bvti(nCv,sdc_nodes), nd_bvpo(nCv,sdc_nodes))
            end if
            write(*,'(a,i2,a,i3,a,es10.2)') 'b2step: SDC active, Radau IIA M = ', &
                sdc_nodes, ', sweeps K = ', sdc_sweeps, ', tol = ', sdc_tol
            write(*,'(a,4f12.8)') '  nodes c    = ', c(1:sdc_nodes)
            write(*,'(a,4f12.8)') '  spacings   = ', dtau(1:sdc_nodes)
            if (sdc_qdelta == 1) write(*,'(a,4f12.8)') &
                '  Q_delta dg = ', (qd(m,m), m = 1, sdc_nodes)
        end if
        call createB2PlasmaSnapshot(nCv, nFc, ns, work)
        call createB2PlasmaSnapshot(nCv, nFc, ns, psnl_entry)
        initialised = .true.
    end subroutine sdc_init

!-----------------------------------------------------------------------
    subroutine reset_ts_factor()
!       b2mndt derives ts_factor from the switches it sees on its first
!       call; internal calls run with b2mndt_ckn = 2, so restore the value
!       implied by the input switches.
        if (input_ckn == 2) then
            ts_factor = 2.0_R8
        else if (input_bdf == 2) then
            ts_factor = 1.5_R8
        else
            ts_factor = 1.0_R8
        end if
    end subroutine reset_ts_factor

!-----------------------------------------------------------------------
    subroutine eval_residual(nout, nCv, nFc, nVx, ns, ismain, ismain0, nscx, &
                             nscxmax, iscx, itim, dtim, ntim, switch, geo, mpg, &
                             st, st_ext, st_avg, snap, respo)
!       F(U) for the current pl: one no-update sweep with the CN residual
!       storage enabled. On exit snap holds the entering state together
!       with its res*0 arrays; pl/dv are restored to the entering state.
        integer :: nout(0:*)
        integer, intent(in) :: nCv, nFc, nVx, ns, ismain, ismain0, nscx, &
                               nscxmax, itim, ntim
        integer :: iscx(0:*)
        real(R8), intent(in) :: dtim
        type(switches), intent(inout) :: switch
        type(geometry), intent(in) :: geo
        type(mapping), intent(inout) :: mpg
        type(B2State), intent(inout) :: st
        type(B2StateExt), intent(inout) :: st_ext
        type(B2Average), intent(inout) :: st_avg
        type(B2PlasmaSnapshot), intent(inout) :: snap
        real(R8), intent(out), optional :: respo(nCv)

        integer :: sv_dummy, sv_nstg(0:2), sv_ckn, sv_moqtlv, sv_moitlv
        integer :: sv_ncall_mndt, sv_ncall_news, sv_itcnt, ierr
        real(R8) :: sv_rxf

        sv_dummy = switch%b2mndt_dummy
        sv_rxf = switch%b2mndt_rxf
        sv_nstg = switch%nstg
        sv_ckn = switch%b2mndt_ckn
        sv_moqtlv = switch%b2mndt_moqtlv
        sv_moitlv = switch%b2mndt_moitlv
        sv_ncall_mndt = ncall_b2mndt
        sv_ncall_news = ncall_b2news_
        sv_itcnt = b2mndt_itcnt

        switch%b2mndt_dummy = 2
        switch%b2mndt_rxf = 0.0_R8
        switch%nstg = 1
        switch%b2mndt_ckn = 2
        switch%b2mndt_moqtlv = -1
        switch%b2mndt_moitlv = -1

!       snapshot of the state to be evaluated; psnl := same state so that the
!       time terms and the b2nxdp potential shift vanish (psnl is restored
!       afterwards: the caller's step measures its time term from it)
        call getB2PlasmaSnapshot(st%pl, st%dv, work)
        psnl_eval_save = st%psnl
        st%psnl = work
        st%psnl%resco0 = 0.0_R8
        st%psnl%resmo0 = 0.0_R8
        st%psnl%reshe0 = 0.0_R8
        st%psnl%reshi0 = 0.0_R8
        st%psnl%reshn0 = 0.0_R8
        st%psnl%reskt0 = 0.0_R8
        st%psnl%reszt0 = 0.0_R8

        call b2mndt(nout, nCv, nFc, nVx, ns, ismain, ismain0, nscx, nscxmax, &
                    iscx, itim, dtim, ntim, switch, geo, mpg, st, st_ext, &
                    st_avg, ierr)

        snap = work
        snap%resco0 = st%dv%resco0
        snap%resmo0 = st%dv%resmo0
        snap%reshe0 = st%dv%reshe0
        snap%reshi0 = st%dv%reshi0
        snap%reshn0 = st%dv%reshn0
        snap%reskt0 = st%dv%reskt0
        snap%reszt0 = st%dv%reszt0
        if (present(respo)) respo = st%dv%respo
        call putB2PlasmaSnapshot(st%pl, st%dv, work)
        st%psnl = psnl_eval_save

        switch%b2mndt_dummy = sv_dummy
        switch%b2mndt_rxf = sv_rxf
        switch%nstg = sv_nstg
        switch%b2mndt_ckn = sv_ckn
        switch%b2mndt_moqtlv = sv_moqtlv
        switch%b2mndt_moitlv = sv_moitlv
        b2mndt_itcnt = sv_itcnt
        ncall_b2mndt = sv_ncall_mndt
        ncall_b2news_ = sv_ncall_news
!       the call may have run the first-call block of b2mndt with ckn = 2
        call reset_ts_factor()
        n_eval_total = n_eval_total + 1
    end subroutine eval_residual

!-----------------------------------------------------------------------
    subroutine mms_update_sources(t, nout, nCv, nFc, nVx, ns, ismain, ismain0, &
                                  nscx, nscxmax, iscx, itim, dtim, ntim, switch, &
                                  geo, mpg, st, st_ext, st_avg)
!       art_* and bv_* := manufactured sources and boundary values at time t
!       by the discrete-operator method. pl/dv are restored on exit.
        real(R8), intent(in) :: t
        integer :: nout(0:*)
        integer, intent(in) :: nCv, nFc, nVx, ns, ismain, ismain0, nscx, &
                               nscxmax, itim, ntim
        integer :: iscx(0:*)
        real(R8), intent(in) :: dtim
        type(switches), intent(inout) :: switch
        type(geometry), intent(in) :: geo
        type(mapping), intent(inout) :: mpg
        type(B2State), intent(inout) :: st
        type(B2StateExt), intent(inout) :: st_ext
        type(B2Average), intent(inout) :: st_avg
        type(B2PlasmaSnapshot), save :: saved, exact
        real(R8) :: respo(nCv)

        call createB2PlasmaSnapshot(nCv, nFc, ns, saved)
        call createB2PlasmaSnapshot(nCv, nFc, ns, exact)
        call getB2PlasmaSnapshot(st%pl, st%dv, saved)
        call mms_set_plasma_state(t, nCv, ns, geo, st%rt%rza, st_ext%ne, st%pl, st%dv)
        art_sna = 0.0_R8
        art_smo = 0.0_R8
        art_she = 0.0_R8
        art_shi = 0.0_R8
        art_shn = 0.0_R8
        art_sch = 0.0_R8
        call mms_set_boundary_values(t, nCv, ns, geo)
        call eval_residual(nout, nCv, nFc, nVx, ns, ismain, ismain0, nscx, &
                           nscxmax, iscx, itim, dtim, ntim, switch, geo, mpg, &
                           st, st_ext, st_avg, exact, respo)
        call mms_fill_sources(t, nCv, ns, geo, mpg, st%rt%rza, st_ext%ne, &
                              exact%resco0, exact%resmo0, exact%reshe0, &
                              exact%reshi0, respo)
        call putB2PlasmaSnapshot(st%pl, st%dv, saved)
        if (sdc_debug_stop /= 0) call xerrab('b2mndt_sdc_debug_stop: stop after first MMS evaluation')
    end subroutine mms_update_sources

!-----------------------------------------------------------------------
    subroutine load_node_sources(m)
        integer, intent(in) :: m
        art_sna = nd_sna(:,:,:,m)
        art_smo = nd_smo(:,:,:,m)
        art_she = nd_she(:,:,m)
        art_shi = nd_shi(:,:,m)
        art_sch = nd_sch(:,:,m)
        art_shn = 0.0_R8
        bv_na = nd_bvna(:,:,m)
        bv_ua = nd_bvua(:,:,m)
        bv_te = nd_bvte(:,m)
        bv_ti = nd_bvti(:,m)
        bv_tn = nd_bvti(:,m)
        bv_po = nd_bvpo(:,m)
    end subroutine load_node_sources

!-----------------------------------------------------------------------
    subroutine node_solve(m, k, nout, nCv, nFc, nVx, ns, ismain, ismain0, nscx, &
                          nscxmax, iscx, itim, dtim, ntim, switch, geo, mpg, &
                          st, st_ext, st_avg, ierr)
!       Solve node m of sweep k with the explicit SDC source b_m in
!       psnl%res*0: backward Euler from node(m-1) over dtau(m)*dtim
!       (sdc_qdelta = 0) or from node(0) over qd(m,m)*dtim (sdc_qdelta = 1).
        integer, intent(in) :: m, k
        integer :: nout(0:*)
        integer, intent(in) :: nCv, nFc, nVx, ns, ismain, ismain0, nscx, &
                               nscxmax, itim, ntim
        integer :: iscx(0:*)
        real(R8), intent(in) :: dtim
        type(switches), intent(inout) :: switch
        type(geometry), intent(in) :: geo
        type(mapping), intent(inout) :: mpg
        type(B2State), intent(inout) :: st
        type(B2StateExt), intent(inout) :: st_ext
        type(B2Average), intent(inout) :: st_avg
        integer, intent(out) :: ierr

        integer :: sv_ckn, sv_nstg2, j, nCi
        real(R8) :: dt_node, w
        logical :: use_ckn

        nCi = mpg%nCi
        if (sdc_qdelta == 0) then
            st%psnl = node(m-1)
        else
            st%psnl = node(0)
        end if
        st%psnl%resco0 = 0.0_R8
        st%psnl%resmo0 = 0.0_R8
        st%psnl%reshe0 = 0.0_R8
        st%psnl%reshi0 = 0.0_R8
        st%psnl%reshn0 = 0.0_R8
        st%psnl%reskt0 = 0.0_R8
        st%psnl%reszt0 = 0.0_R8
        if (sdc_qdelta == 0) then
            if (k > 0) then
!               b_m = (1/dtau_m) sum_j S_mj F^k_j - F^k_m over interior cells
                do j = 1, sdc_nodes
                    w = s(m,j)/dtau(m)
                    call add_fk(j, w)
                end do
                call add_fk(m, -1.0_R8)
            end if
            use_ckn = k > 0
            dt_node = dtau(m)*dtim
        else
!           b_m = (1/qd_mm) [ sum_j Q_mj F^k_j + sum_{j<m} qd_mj (F^{k+1}_j - F^k_j) ]
!                 - F^k_m over interior cells; F^k = 0 in sweep 0
            if (k > 0) then
                do j = 1, sdc_nodes
                    w = q(m,j)/qd(m,m)
                    call add_fk(j, w)
                end do
                call add_fk(m, -1.0_R8)
            end if
            do j = 1, m - 1
                w = qd(m,j)/qd(m,m)
                call add_node_f(j, w)
                if (k > 0) call add_fk(j, -w)
            end do
            use_ckn = k > 0 .or. m > 1
            dt_node = qd(m,m)*dtim
        end if

        sv_ckn = switch%b2mndt_ckn
        sv_nstg2 = switch%nstg(2)
!       the CN channel carries b_m; node 1 of sweep 0 has b_m = 0 and runs as
!       plain implicit Euler, so that the first-call block of b2mndt (which
!       derives ts_factor from b2mndt_ckn) sees the input value
        if (use_ckn) switch%b2mndt_ckn = 2
        if (sdc_nstg2 > 0) switch%nstg(2) = sdc_nstg2
        call b2mndt(nout, nCv, nFc, nVx, ns, ismain, ismain0, nscx, nscxmax, &
                    iscx, itim, dt_node, ntim, switch, geo, mpg, st, st_ext, &
                    st_avg, ierr)
        switch%b2mndt_ckn = sv_ckn
        switch%nstg(2) = sv_nstg2

        if (ierr == 0) then
            if (.not. (all(ieee_is_finite(st%pl%na(1:nCi,:))) .and. &
                       all(ieee_is_finite(st%pl%ua(1:nCi,:))) .and. &
                       all(ieee_is_finite(st%pl%te(1:nCi))) .and. &
                       all(ieee_is_finite(st%pl%ti(1:nCi))) .and. &
                       all(ieee_is_finite(st%pl%po(1:nCi))) .and. &
                       all(st%pl%na(1:nCi,:) > 0.0_R8) .and. &
                       all(st%pl%te(1:nCi) > 0.0_R8) .and. &
                       all(st%pl%ti(1:nCi) > 0.0_R8))) ierr = 1
        end if
        if (ierr == 0) call getB2PlasmaSnapshot(st%pl, st%dv, node(m))

    contains

        subroutine add_fk(j, w)
!           psnl%res*0 += w F^k_j (previous sweep) over interior cells
            integer, intent(in) :: j
            real(R8), intent(in) :: w
            st%psnl%resco0(1:nCi,:) = st%psnl%resco0(1:nCi,:) + w*fk_co(1:nCi,:,j)
            st%psnl%resmo0(1:nCi,:) = st%psnl%resmo0(1:nCi,:) + w*fk_mo(1:nCi,:,j)
            st%psnl%reshe0(1:nCi) = st%psnl%reshe0(1:nCi) + w*fk_he(1:nCi,j)
            st%psnl%reshi0(1:nCi) = st%psnl%reshi0(1:nCi) + w*fk_hi(1:nCi,j)
        end subroutine add_fk

        subroutine add_node_f(j, w)
!           psnl%res*0 += w F^{k+1}_j (current sweep, node(j)%res*0)
            integer, intent(in) :: j
            real(R8), intent(in) :: w
            st%psnl%resco0(1:nCi,:) = st%psnl%resco0(1:nCi,:) + w*node(j)%resco0(1:nCi,:)
            st%psnl%resmo0(1:nCi,:) = st%psnl%resmo0(1:nCi,:) + w*node(j)%resmo0(1:nCi,:)
            st%psnl%reshe0(1:nCi) = st%psnl%reshe0(1:nCi) + w*node(j)%reshe0(1:nCi)
            st%psnl%reshi0(1:nCi) = st%psnl%reshi0(1:nCi) + w*node(j)%reshi0(1:nCi)
        end subroutine add_node_f

    end subroutine node_solve

!-----------------------------------------------------------------------
    function sdc_residual(ns, nCi, dtim, geo) result(res)
!       Collocation residual at the last node,
!       R = U_M - U_0 - dt sum_j Q_Mj F_j (volume-integrated conserved
!       variables), as an L2 norm relative to the step increment U_M - U_0,
!       maximum over the equations.
        integer, intent(in) :: ns, nCi
        real(R8), intent(in) :: dtim
        type(geometry), intent(in) :: geo
        real(R8) :: res
        real(R8) :: u0(nCi), um(nCi), r(nCi), f(nCi)
        integer :: is, j, mm
        real(R8) :: rn, dn

        mm = sdc_nodes
        res = 0.0_R8
        do is = 0, ns - 1
!           particles
            u0 = geo%cvVol(1:nCi)*node(0)%na(1:nCi,is)
            um = geo%cvVol(1:nCi)*node(mm)%na(1:nCi,is)
            f = 0.0_R8
            do j = 1, mm
                f = f + q(mm,j)*node(j)%resco0(1:nCi,is)
            end do
            call accumulate()
!           momentum (am*mp cancels in the ratio only if included in both; use conserved form)
            u0 = geo%cvVol(1:nCi)*geo%cvHz(1:nCi)*mass(is)* &
                 node(0)%na(1:nCi,is)*node(0)%ua(1:nCi,is)
            um = geo%cvVol(1:nCi)*geo%cvHz(1:nCi)*mass(is)* &
                 node(mm)%na(1:nCi,is)*node(mm)%ua(1:nCi,is)
            f = 0.0_R8
            do j = 1, mm
                f = f + q(mm,j)*node(j)%resmo0(1:nCi,is)
            end do
            call accumulate()
        end do
!       electron energy
        u0 = 1.5_R8*geo%cvVol(1:nCi)*node(0)%ne(1:nCi)*node(0)%te(1:nCi)
        um = 1.5_R8*geo%cvVol(1:nCi)*node(mm)%ne(1:nCi)*node(mm)%te(1:nCi)
        f = 0.0_R8
        do j = 1, mm
            f = f + q(mm,j)*node(j)%reshe0(1:nCi)
        end do
        call accumulate()
!       ion (+ neutral, tn_style = 0) energy
        u0 = 1.5_R8*geo%cvVol(1:nCi)*node(0)%ni(1:nCi,0)*node(0)%ti(1:nCi)
        um = 1.5_R8*geo%cvVol(1:nCi)*node(mm)%ni(1:nCi,0)*node(mm)%ti(1:nCi)
        f = 0.0_R8
        do j = 1, mm
            f = f + q(mm,j)*node(j)%reshi0(1:nCi)
        end do
        call accumulate()

    contains

        subroutine accumulate()
            r = um - u0 - dtim*f
            rn = sqrt(sum(r**2))
            dn = sqrt(sum((um - u0)**2))
            if (dn > 0.0_R8) then
                res = max(res, rn/dn)
            else if (rn > 0.0_R8) then
                res = max(res, huge(1.0_R8))
            end if
        end subroutine accumulate

    end function sdc_residual

!-----------------------------------------------------------------------
    function mass(is) result(m)
        use b2mod_constants, only: mp
        use b2mod_b2cmpa, only: am
        integer, intent(in) :: is
        real(R8) :: m
        m = am(is)*mp
    end function mass

end module b2mod_sdc
