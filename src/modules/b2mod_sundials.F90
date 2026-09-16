module b2mod_sundials
    use b2mod_types, only: R8
    use b2mod_switches, only: switches
    use b2mod_ad, only: ncall_b2news_
    use b2mod_numerics_namelist, only: min_na
    use b2us_geo, only: geometry
    use b2us_map, only: mapping
    use b2us_plasma, only: B2State, B2StateExt, B2Average
#if defined(USE_SUNDIALS)
    use, intrinsic :: iso_c_binding, only: &
        c_null_ptr, c_ptr, c_int, c_int64_t, c_long, c_double, &
        c_associated, c_funloc
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fsundials_core_mod, only: & ! IGNORE
        SUN_COMM_NULL, &
        N_Vector, &
        FSUNContext_Free, FSUNContext_Create, &
        FN_VDestroy, FN_VConst, FN_VGetArrayPointer
    use fnvector_serial_mod, only: & ! IGNORE
        FN_VNew_Serial
    use fkinsol_mod, only: & ! IGNORE
        KIN_FP, KIN_SUCCESS, KIN_MAXITER_REACHED, &
        FKINSol, &
        FKINFree, FKINCreate, FKINInit, FKINSetNumMaxIters, FKINSetMAA, &
        FKINSetDamping, FKINSetDampingAA, &
        FKINSetFuncNormTol, FKINSetScaledStepTol, FKINGetNumNonlinSolvIters
#endif

    implicit none
    public

!   Status returned to b2mndt. The outcome of a solve is decided by the
!   wrapper's own residual bookkeeping (see residual_is_acceptable), not by
!   the KINSOL return flag, which is reported for information only.
    integer, parameter :: B2_KINSOL_FATAL = -1     ! re-entrancy or SUNDIALS bookkeeping failure: abort
    integer, parameter :: B2_KINSOL_SUCCESS = 0    ! residual target reached
    integer, parameter :: B2_KINSOL_FALLBACK = 1   ! nothing acceptable: entry state restored, caller sweeps
    integer, parameter :: B2_KINSOL_MAXITER = 2    ! budget exhausted, best iterate accepted by the accept factor
    integer, parameter :: B2_KINSOL_PARTIAL = 3    ! KINSOL stopped early, best iterate accepted by the accept factor

!   Per-field residual slots (transformed L-infinity over interior cells).
    integer, parameter :: B2_KINSOL_NFIELD = 8
    integer, parameter :: B2_KINSOL_FIELD_NA = 1, B2_KINSOL_FIELD_UA = 2, &
        B2_KINSOL_FIELD_PO = 3, B2_KINSOL_FIELD_TE = 4, &
        B2_KINSOL_FIELD_TI = 5, B2_KINSOL_FIELD_TN = 6, &
        B2_KINSOL_FIELD_KT = 7, B2_KINSOL_FIELD_ZT = 8
    character(len=2), parameter :: B2_KINSOL_FIELD_NAME(B2_KINSOL_NFIELD) = &
        ['na', 'ua', 'po', 'te', 'ti', 'tn', 'kt', 'zt']

#if defined(USE_SUNDIALS)
    logical, save :: solver_active = .false.
    logical, save :: active_include_po = .true.
    logical, save :: active_include_tn = .true.
    logical, save :: active_include_kt = .false.
    logical, save :: active_include_zt = .false.
    logical, save :: best_output_valid = .false.
    logical, save :: target_reached = .false.
    integer, save :: n_active_callback = 0
    integer, save :: n_active_sync = 0
    integer, save :: n_active_projection = 0
    integer, save :: best_evaluation = 0
    integer, save :: best_worst_index = 0
    integer, save :: active_iout = 0
    integer, save :: active_ncall_b2news = 0
    integer, save :: active_first_species = 0
    integer, save :: active_last_species = -1
    integer, save :: active_nCv = 0
    integer, save :: active_nCi = 0
    integer, save :: active_nFc = 0
    integer, save :: active_nVx = 0
    integer, save :: active_ns = 0
    integer, save :: active_nscx = 0
    integer, save :: active_nscxmax = 0
    integer, save :: active_ismain = 0
    integer, save :: active_ismain0 = 0
    integer, save :: active_state_size = 0
    integer, save :: active_sral = 0
    integer, save :: workspace_state_size = 0
    integer, allocatable, save :: active_iscx(:)
    real(R8), save :: active_dtim = 0.0_R8
    real(R8), save :: active_trust_radius = 0.0_R8
    real(R8), save :: active_rxf = 0.0_R8
    real(R8), save :: active_fnorm_floor = 0.0_R8
    real(R8), save :: active_fnorm_rtol = 0.0_R8
    real(R8), save :: active_fnorm_target = 0.0_R8
    real(R8), save :: active_accept_factor = 0.0_R8
    real(R8), save :: active_residual_first = 0.0_R8
    real(R8), save :: active_residual_min = 0.0_R8
    real(R8), save :: active_residual_last = 0.0_R8
    real(R8), save :: best_field_residual(B2_KINSOL_NFIELD) = 0.0_R8
    real(R8), allocatable, save :: scale_na(:), scale_ua(:)
    real(R8), allocatable, save :: unpack_buffer(:)
    real(R8), save :: scale_po = 1.0_R8
    real(R8), save :: scale_te = 1.0_R8
    real(R8), save :: scale_ti = 1.0_R8
    real(R8), save :: scale_tn = 1.0_R8
    real(R8), save :: scale_kt = 1.0_R8
    real(R8), save :: scale_zt = 1.0_R8
!   interior_mask(k) is true for state-vector entries that belong to an
!   interior cell (1..nCi). Guard cells stay in the vector -- they are inputs
!   to b2news_m and are restored exactly on fallback -- but they are excluded
!   from every norm the solver steers by.
    logical, allocatable, save :: interior_mask(:)
    real(c_double), allocatable, save :: initial_state_value(:)
    real(c_double), allocatable, save :: best_output_value(:)
    real(c_double), allocatable, save :: projected_input_value(:)
    type(switches), pointer, save :: active_switch => null()
    type(geometry), pointer, save :: active_geo => null()
    type(mapping), pointer, save :: active_mpg => null()
    type(B2State), pointer, save :: active_st => null()
    type(B2StateExt), pointer, save :: active_st_ext => null()
    type(B2Average), pointer, save :: active_st_avg => null()
    type(c_ptr), save :: sundials_context = c_null_ptr
    type(c_ptr), save :: kinsol_memory = c_null_ptr
    type(N_Vector), pointer, save :: workspace_f_scale_vector => null()
    type(N_Vector), pointer, save :: workspace_u_scale_vector => null()
    type(N_Vector), pointer, save :: workspace_state_vector => null()
#endif

contains

    subroutine b2_sundials_solve(nCv, nFc, nVx, ns, nscx, &
        iscx, nscxmax, ismain, ismain0, dtim, switch, geo, mpg, st, &
        st_ext, st_avg, max_iteration, maa, fnorm_tolerance, &
        fnorm_rtolerance, step_tolerance, &
        damping, trust_radius, iteration_rxf, &
        iteration_sral, accept_factor, damping_aa, iteration_iout, &
        n_iteration, n_evaluation, n_sync_evaluation, n_projection, &
        n_best_evaluation, &
        residual_first, residual_min, residual_last, field_residual, &
        worst_field, worst_species, worst_cell, &
        status, kinsol_flag)
        integer, intent(in) :: nCv, nFc, nVx, ns, nscx, nscxmax
        integer, intent(in) :: iscx(0:nscxmax-1), ismain, ismain0
        real(R8), intent(in) :: dtim
        type(switches), target, intent(inout) :: switch
        type(geometry), target, intent(in) :: geo
        type(mapping), target, intent(inout) :: mpg
        type(B2State), target, intent(inout) :: st
        type(B2StateExt), target, intent(inout) :: st_ext
        type(B2Average), target, intent(in) :: st_avg
        integer, intent(in) :: max_iteration, maa
!       fnorm_tolerance: absolute floor of the residual target.
!       fnorm_rtolerance: target = max(floor, rtol * res_first); 0 = floor only.
!       step_tolerance: passed to FKINSetScaledStepTol for completeness only.
!         KIN_FP in the linked SUNDIALS tests mxiter and fnormtol only, so
!         this value has no effect on the fixed-point iteration.
!       accept_factor: a best iterate with res_min <= accept_factor * res_first
!         is accepted even if the target was not reached; 0 = target only.
!       damping: damping of the plain fixed-point step (KINSetDamping).
!       damping_aa: damping of the Anderson-accelerated step
!         (KINSetDampingAA); <= 0 means the same value as damping.
!       iteration_iout: 1 prints one line per map evaluation.
!       worst_*: field slot, species (-1 if not a species field) and cell of
!         the state-vector entry with the largest residual on the best
!         evaluation.
        real(R8), intent(in) :: fnorm_tolerance, fnorm_rtolerance, &
            step_tolerance, damping, trust_radius, iteration_rxf, &
            accept_factor, damping_aa
        integer, intent(in) :: iteration_sral, iteration_iout
        integer, intent(out) :: n_iteration, n_evaluation, n_sync_evaluation
        integer, intent(out) :: n_projection, n_best_evaluation
        real(R8), intent(out) :: residual_first, residual_min, residual_last
        real(R8), intent(out) :: field_residual(B2_KINSOL_NFIELD)
        integer, intent(out) :: worst_field, worst_species, worst_cell
        integer, intent(out) :: status, kinsol_flag

#if defined(USE_SUNDIALS)
        integer(c_int) :: flag, cleanup_flag
        integer(c_long) :: n_nonlinear_iteration(1), maa_value
        real(c_double), pointer :: state_value(:)
        real(R8) :: effective_fnorm_tolerance, effective_step_tolerance
        real(R8) :: effective_damping_aa
        logical :: best_is_current, valid, accepted

        n_iteration = 0
        n_evaluation = 0
        n_sync_evaluation = 0
        n_projection = 0
        n_best_evaluation = 0
        residual_first = 0.0_R8
        residual_min = 0.0_R8
        residual_last = 0.0_R8
        field_residual = 0.0_R8
        worst_field = 0
        worst_species = -1
        worst_cell = 0
        best_worst_index = 0
        active_iout = iteration_iout
        active_residual_first = 0.0_R8
        active_residual_min = 0.0_R8
        active_residual_last = 0.0_R8
        active_fnorm_target = 0.0_R8
        best_field_residual = 0.0_R8
        active_trust_radius = trust_radius
        active_rxf = iteration_rxf
        active_sral = iteration_sral
        active_fnorm_rtol = fnorm_rtolerance
        active_accept_factor = accept_factor
        status = B2_KINSOL_FALLBACK
        kinsol_flag = 0

        if (solver_active) then
            status = B2_KINSOL_FATAL
            kinsol_flag = -1
            return
        end if
        if (max_iteration .lt. 1) return

        call configure_active_state(switch, mpg, st, valid)
        if (.not. valid) goto 900
!       Put na inside the window b2news_m clamps it to before anything
!       measures the state: otherwise the scales and res_first include the
!       distance back to [na_min, na_max], which is unbounded in log
!       variables and is not a property of the map.
        call clamp_state_na(switch, mpg, st)
        call initialize_scale(st, valid)
        if (.not. valid) goto 900
        call ensure_active_storage(nscxmax, valid)
        if (.not. valid) goto 900
        active_iscx = iscx

        call ensure_workspace(active_state_size, kinsol_flag, valid)
        if (.not. valid) goto 900
        state_value => FN_VGetArrayPointer(workspace_state_vector)
        if (.not. associated(state_value)) then
            kinsol_flag = -1
            goto 900
        end if
        call pack_state(st, state_value, valid)
        if (.not. valid) goto 900
        initial_state_value = state_value
        call FN_VConst(1.0_c_double, workspace_u_scale_vector)
        call FN_VConst(1.0_c_double, workspace_f_scale_vector)

        flag = FKINInit(kinsol_memory, &
            c_funloc(b2_sundials_rhs_fixed_point), workspace_state_vector)
        if (flag .ne. KIN_SUCCESS) then
            kinsol_flag = flag
            goto 900
        end if

        maa_value = int(min(maa, max(0, max_iteration - 1)), c_long)
        flag = FKINSetNumMaxIters(kinsol_memory, int(max_iteration, c_long))
        if (flag .ne. KIN_SUCCESS) then
            kinsol_flag = flag
            goto 900
        end if
        flag = FKINSetMAA(kinsol_memory, maa_value)
        if (flag .ne. KIN_SUCCESS) then
            kinsol_flag = flag
            goto 900
        end if
        flag = FKINSetDamping(kinsol_memory, real(damping, c_double))
        if (flag .ne. KIN_SUCCESS) then
            kinsol_flag = flag
            goto 900
        end if
        effective_damping_aa = damping_aa
        if (effective_damping_aa .le. 0.0_R8) effective_damping_aa = damping
        flag = FKINSetDampingAA(kinsol_memory, &
            real(effective_damping_aa, c_double))
        if (flag .ne. KIN_SUCCESS) then
            kinsol_flag = flag
            goto 900
        end if
!       KINSOL's own KIN_FP convergence test measures the *update*
!       ||u_new - u||, which is damping * ||G(u) - u|| for a plain
!       fixed-point step and can be arbitrarily small for a degenerate
!       Anderson step while the residual is not. It is therefore disabled
!       here (a positive tolerance far below anything reachable; 0 would
!       select the SUNDIALS default) and the callback stops the iteration
!       itself when the true residual ||G(u) - u|| reaches the target.
        flag = FKINSetFuncNormTol(kinsol_memory, tiny(1.0_c_double))
        if (flag .ne. KIN_SUCCESS) then
            kinsol_flag = flag
            goto 900
        end if
        effective_fnorm_tolerance = fnorm_tolerance
        if (effective_fnorm_tolerance .le. 0.0_R8) then
            effective_fnorm_tolerance = epsilon(1.0_R8)**(1.0_R8/3.0_R8)
        end if
        active_fnorm_floor = effective_fnorm_tolerance
        active_fnorm_target = effective_fnorm_tolerance
!       No effect in KIN_FP (see the argument note); set for completeness.
        effective_step_tolerance = step_tolerance
        if (effective_step_tolerance .le. 0.0_R8) then
            effective_step_tolerance = epsilon(1.0_R8)**(2.0_R8/3.0_R8)
        end if
        flag = FKINSetScaledStepTol(kinsol_memory, &
            real(effective_step_tolerance, c_double))
        if (flag .ne. KIN_SUCCESS) then
            kinsol_flag = flag
            goto 900
        end if

        active_nFc = nFc
        active_nVx = nVx
        active_ns = ns
        active_nscx = nscx
        active_nscxmax = nscxmax
        active_ismain = ismain
        active_ismain0 = ismain0
        active_dtim = dtim
        n_active_callback = 0
        n_active_sync = 0
        n_active_projection = 0
        best_evaluation = 0
        best_output_valid = .false.
        target_reached = .false.
        active_ncall_b2news = ncall_b2news_
        active_switch => switch
        active_geo => geo
        active_mpg => mpg
        active_st => st
        active_st_ext => st_ext
        active_st_avg => st_avg
        solver_active = .true.

        flag = FKINSol(kinsol_memory, workspace_state_vector, &
            KIN_FP, workspace_u_scale_vector, workspace_f_scale_vector)
        kinsol_flag = flag
        n_evaluation = n_active_callback
        cleanup_flag = FKINGetNumNonlinSolvIters(kinsol_memory, &
            n_nonlinear_iteration)
        if (cleanup_flag == KIN_SUCCESS) then
            n_iteration = int(n_nonlinear_iteration(1))
        else if (n_active_callback .gt. 0) then
            status = B2_KINSOL_FATAL
            kinsol_flag = cleanup_flag
            goto 900
        else
            goto 900
        end if

!       The state handed back is always the output of the evaluation with the
!       smallest residual, i.e. G(u*) for the u* that came closest to a fixed
!       point -- a state that a real sweep produced. KINSOL's own returned
!       vector is not used: it is an accelerated combination whose residual
!       has never been evaluated, and it may be a KIN_SYSFUNC_FAIL stop that
!       the callback requested on purpose (target reached).
        accepted = residual_is_acceptable()
        if (.not. accepted) then
            status = B2_KINSOL_FALLBACK
            goto 900
        end if

        best_is_current = best_evaluation == n_active_callback
        if (.not. best_is_current) then
            call unpack_state(best_output_value, st, valid)
            if (.not. valid) then
                status = B2_KINSOL_FALLBACK
                goto 900
            end if
            call sync_derived_state(nCv, nFc, nVx, ns, nscx, iscx, nscxmax, &
                ismain, ismain0, dtim, switch, geo, mpg, st, st_ext, st_avg, &
                iteration_sral)
        end if

        if (target_reached) then
            status = B2_KINSOL_SUCCESS
        else if (flag == KIN_MAXITER_REACHED) then
            status = B2_KINSOL_MAXITER
        else
            status = B2_KINSOL_PARTIAL
        end if

900 continue
        if (solver_active .and. &
            (status == B2_KINSOL_FATAL .or. &
                status == B2_KINSOL_FALLBACK)) then
!           Restore the entry state. The vector holds every cell, guard cells
!           included, so the plasma fields come back exactly up to the
!           log/exp round trip; derived quantities are rebuilt by the sync.
            call unpack_state(initial_state_value, st, valid)
            if (valid .and. n_active_callback .gt. 0) then
                call sync_derived_state(nCv, nFc, nVx, ns, nscx, iscx, &
                    nscxmax, ismain, ismain0, dtim, switch, geo, mpg, st, &
                    st_ext, st_avg, iteration_sral)
            end if
        end if
        residual_first = active_residual_first
        residual_min = active_residual_min
        residual_last = active_residual_last
        field_residual = best_field_residual
        n_best_evaluation = best_evaluation
        if (best_worst_index .gt. 0) call decode_index(best_worst_index, &
            worst_field, worst_species, worst_cell)
        if (n_active_callback .gt. 0) then
            ncall_b2news_ = active_ncall_b2news + n_active_callback + &
                n_active_sync
        end if
        n_evaluation = n_active_callback
        n_sync_evaluation = n_active_sync
        n_projection = n_active_projection
        call clear_active_state()
#else
        n_iteration = 0
        n_evaluation = 0
        n_sync_evaluation = 0
        n_projection = 0
        n_best_evaluation = 0
        residual_first = 0.0_R8
        residual_min = 0.0_R8
        residual_last = 0.0_R8
        field_residual = 0.0_R8
        worst_field = 0
        worst_species = -1
        worst_cell = 0
        status = B2_KINSOL_FATAL
        kinsol_flag = 0
#endif
    end subroutine b2_sundials_solve

#if defined(USE_SUNDIALS)
!   One b2news_m call with rxf = 0: every block solve is set up and every
!   derived quantity (fluxes, ne/ni/nn, transport coefficients, residuals,
!   sources) is rebuilt for the current plasma state, but no correction is
!   applied, so the state does not move.
    subroutine sync_derived_state(nCv, nFc, nVx, ns, nscx, iscx, nscxmax, &
        ismain, ismain0, dtim, switch, geo, mpg, st, st_ext, st_avg, &
        iteration_sral)
        integer, intent(in) :: nCv, nFc, nVx, ns, nscx, nscxmax
        integer, intent(in) :: iscx(0:nscxmax-1), ismain, ismain0
        real(R8), intent(in) :: dtim
        type(switches), intent(inout) :: switch
        type(geometry), intent(in) :: geo
        type(mapping), intent(inout) :: mpg
        type(B2State), intent(inout) :: st
        type(B2StateExt), intent(inout) :: st_ext
        type(B2Average), intent(in) :: st_avg
        integer, intent(in) :: iteration_sral
        real(R8) :: saved_rxf
        integer :: saved_sral

        saved_rxf = switch%b2mndt_rxf
        saved_sral = switch%no_b2sral_call
        switch%b2mndt_rxf = 0.0_R8
        if (iteration_sral .eq. 1) switch%no_b2sral_call = 0
        ncall_b2news_ = active_ncall_b2news
        call b2news_m(nCv, nFc, nVx, ns, nscx, iscx, nscxmax, ismain, &
            ismain0, dtim, switch, geo, mpg, st, st_ext, st_avg, .false.)
        switch%b2mndt_rxf = saved_rxf
        switch%no_b2sral_call = saved_sral
        n_active_sync = n_active_sync + 1
    end subroutine sync_derived_state
#endif

#if defined(USE_SUNDIALS)
!   The best evaluated iterate is acceptable if its residual reached the
!   target, or fell to accept_factor times the first residual. An iterate
!   that never improved on res_first is never acceptable: for that one the
!   caller's ordinary sweep is the better use of the step.
    logical function residual_is_acceptable() result(acceptable)
        acceptable = .false.
        if (.not. best_output_valid) return
        if (n_active_callback .lt. 1) return
        if (target_reached) then
            acceptable = .true.
            return
        end if
        if (active_accept_factor .gt. 0.0_R8 .and. &
            active_residual_first .gt. 0.0_R8) then
            acceptable = active_residual_min .le. &
                active_accept_factor * active_residual_first .and. &
                active_residual_min .lt. active_residual_first
        end if
    end function residual_is_acceptable
#endif

#if defined(USE_SUNDIALS)
    integer(c_int) function b2_sundials_rhs_fixed_point(sunvec_in, &
        sunvec_out, user_data) result(flag) bind(C)
        type(N_Vector) :: sunvec_in, sunvec_out
        type(c_ptr), value :: user_data
        real(c_double), pointer :: value_in(:), value_out(:)
        real(R8) :: change, saved_rxf
        real(R8) :: field_change(B2_KINSOL_NFIELD)
        integer :: saved_sral, worst_index, field, species, cell
        logical :: valid, projected

        flag = -1_c_int
        if (.not. solver_active) return
        value_in => FN_VGetArrayPointer(sunvec_in)
        value_out => FN_VGetArrayPointer(sunvec_out)
        if (.not. associated(value_in) .or. &
            .not. associated(value_out)) return
        if (size(value_in) .ne. active_state_size .or. &
            size(value_out) .ne. active_state_size) return
        if (.not. all(ieee_is_finite(value_in))) return

!       Trust region: instead of failing the step (KIN_FP has no recovery
!       from a callback error, so a rejection would discard every
!       evaluation already spent), project the input onto the L-infinity
!       ball of radius trust around the entry state and evaluate the map
!       there. The iteration then runs on G(P(u)), which has the same fixed
!       points as G inside the ball. The residual is still measured against
!       the unprojected input so that KINSOL's iteration stays consistent.
        projected = .false.
        if (active_trust_radius .gt. 0.0_R8) then
            projected_input_value = min(max(value_in, &
                initial_state_value - real(active_trust_radius, c_double)), &
                initial_state_value + real(active_trust_radius, c_double))
            projected = any(projected_input_value .ne. value_in)
            if (projected) n_active_projection = n_active_projection + 1
            call unpack_state(projected_input_value, active_st, valid)
        else
            call unpack_state(value_in, active_st, valid)
        end if
        if (.not. valid) return

        n_active_callback = n_active_callback + 1
        ncall_b2news_ = active_ncall_b2news
        saved_rxf = active_switch%b2mndt_rxf
        saved_sral = active_switch%no_b2sral_call
        if (active_rxf .gt. 0.0_R8) active_switch%b2mndt_rxf = active_rxf
        if (active_sral .eq. 1) active_switch%no_b2sral_call = 0
        call b2news_m(active_nCv, active_nFc, active_nVx, active_ns, &
            active_nscx, active_iscx, active_nscxmax, active_ismain, &
            active_ismain0, active_dtim, active_switch, active_geo, &
            active_mpg, active_st, active_st_ext, active_st_avg, .false.)
        active_switch%b2mndt_rxf = saved_rxf
        active_switch%no_b2sral_call = saved_sral
        ncall_b2news_ = active_ncall_b2news

        call pack_state(active_st, value_out, valid)
        if (.not. valid) return

!       Residual ||G(u) - u|| in transformed units over interior cells only.
        call masked_residual(value_out, value_in, change, field_change, &
            worst_index)
        active_residual_last = change
        if (active_iout .ne. 0) then
            call decode_index(worst_index, field, species, cell)
            write(*, '(a,i4,a,es9.2,a,l1,3a,i3,a,i6)') &
                ' KINSOL: eval ', n_active_callback, ' res ', change, &
                ' proj ', projected, ' worst ', &
                B2_KINSOL_FIELD_NAME(max(1, field)), ' is ', species, &
                ' cell ', cell
        end if
        if (n_active_callback .eq. 1) then
            active_residual_first = change
            active_fnorm_target = active_fnorm_floor
            if (active_fnorm_rtol .gt. 0.0_R8) then
                active_fnorm_target = max(active_fnorm_floor, &
                    active_fnorm_rtol * change)
            end if
        end if
        if (n_active_callback .eq. 1 .or. change .lt. active_residual_min) then
            active_residual_min = change
            best_output_value = value_out
            best_field_residual = field_change
            best_worst_index = worst_index
            best_evaluation = n_active_callback
            best_output_valid = .true.
        end if

        if (change .le. active_fnorm_target) then
!           Target reached. KIN_FP offers no way for the callback to declare
!           convergence, so stop the iteration with a negative return: KINFP
!           breaks out with KIN_SYSFUNC_FAIL and leaves its vector alone;
!           b2_sundials_solve ignores that vector and uses best_output_value.
            target_reached = .true.
            flag = -1_c_int
            return
        end if

        flag = 0_c_int
    end function b2_sundials_rhs_fixed_point
#endif

#if defined(USE_SUNDIALS)
!   L-infinity norm of (a - b) over interior entries, total and per field.
!   Walks the pack_state layout: na per species, ua per species, po, te, ti,
!   tn, kt, zt, each block active_nCv long.
    subroutine masked_residual(a, b, total, per_field, worst_index)
        real(c_double), intent(in) :: a(:), b(:)
        real(R8), intent(out) :: total, per_field(B2_KINSOL_NFIELD)
        integer, intent(out) :: worst_index
        integer :: first, last, is, field, block_arg
        real(R8) :: block_max

        per_field = 0.0_R8
        total = 0.0_R8
        worst_index = 0
        first = 1
        do is = active_first_species, active_last_species
            call next_block(B2_KINSOL_FIELD_NA)
        end do
        do is = active_first_species, active_last_species
            call next_block(B2_KINSOL_FIELD_UA)
        end do
        if (active_include_po) call next_block(B2_KINSOL_FIELD_PO)
        call next_block(B2_KINSOL_FIELD_TE)
        call next_block(B2_KINSOL_FIELD_TI)
        if (active_include_tn) call next_block(B2_KINSOL_FIELD_TN)
        if (active_include_kt) call next_block(B2_KINSOL_FIELD_KT)
        if (active_include_zt) call next_block(B2_KINSOL_FIELD_ZT)

    contains

        subroutine next_block(slot)
            integer, intent(in) :: slot
            field = slot
            last = first + active_nCv - 1
            block_arg = maxloc(abs(real(a(first:last), R8) - &
                real(b(first:last), R8)), dim=1, &
                mask=interior_mask(first:last))
            if (block_arg .ge. 1) then
                block_max = abs(real(a(first+block_arg-1), R8) - &
                    real(b(first+block_arg-1), R8))
                per_field(field) = max(per_field(field), block_max)
                if (block_max .gt. total .or. worst_index == 0) then
                    total = block_max
                    worst_index = first + block_arg - 1
                end if
            end if
            first = last + 1
        end subroutine next_block
    end subroutine masked_residual
#endif

#if defined(USE_SUNDIALS)
!   Field slot, species (-1 for non-species fields) and cell index of a
!   state-vector entry, following the pack_state block order.
    subroutine decode_index(index, field, species, cell)
        integer, intent(in) :: index
        integer, intent(out) :: field, species, cell
        integer :: iblock, k, n_species

        field = 0
        species = -1
        cell = 0
        if (index .lt. 1 .or. active_nCv .lt. 1) return
        iblock = (index - 1) / active_nCv + 1
        cell = index - (iblock - 1) * active_nCv
        n_species = active_last_species - active_first_species + 1
        if (iblock .le. n_species) then
            field = B2_KINSOL_FIELD_NA
            species = active_first_species + iblock - 1
            return
        end if
        if (iblock .le. 2 * n_species) then
            field = B2_KINSOL_FIELD_UA
            species = active_first_species + iblock - n_species - 1
            return
        end if
        k = 2 * n_species
        if (active_include_po) then
            k = k + 1
            if (iblock == k) then
                field = B2_KINSOL_FIELD_PO
                return
            end if
        end if
        k = k + 1
        if (iblock == k) then
            field = B2_KINSOL_FIELD_TE
            return
        end if
        k = k + 1
        if (iblock == k) then
            field = B2_KINSOL_FIELD_TI
            return
        end if
        if (active_include_tn) then
            k = k + 1
            if (iblock == k) then
                field = B2_KINSOL_FIELD_TN
                return
            end if
        end if
        if (active_include_kt) then
            k = k + 1
            if (iblock == k) then
                field = B2_KINSOL_FIELD_KT
                return
            end if
        end if
        if (active_include_zt) then
            k = k + 1
            if (iblock == k) field = B2_KINSOL_FIELD_ZT
        end if
    end subroutine decode_index
#endif

#if defined(USE_SUNDIALS)
    subroutine clear_active_state()
        solver_active = .false.
        n_active_callback = 0
        n_active_sync = 0
        n_active_projection = 0
        best_evaluation = 0
        best_worst_index = 0
        best_output_valid = .false.
        target_reached = .false.
        nullify(active_switch, active_geo, active_mpg, active_st, &
            active_st_ext, active_st_avg)
    end subroutine clear_active_state
#endif

#if defined(USE_SUNDIALS)
    subroutine configure_active_state(switch, mpg, st, valid)
        type(switches), intent(in) :: switch
        type(mapping), intent(in) :: mpg
        type(B2State), intent(in) :: st
        logical, intent(out) :: valid
        integer :: n_active_species, n_block, first, last, k, &
            allocation_status

        active_first_species = max(lbound(st%pl%na, 2), switch%nsmin)
        active_last_species = min(ubound(st%pl%na, 2), switch%nsmax - 1)
        valid = active_last_species .ge. active_first_species
        if (.not. valid) return
        valid = lbound(st%pl%ua, 2) .le. active_first_species .and. &
            ubound(st%pl%ua, 2) .ge. active_last_species
        if (.not. valid) return

!       Every field block is one full cell array; all blocks must agree.
        active_nCv = size(st%pl%na, 1)
        valid = active_nCv .gt. 0 .and. &
            size(st%pl%ua, 1) == active_nCv .and. &
            size(st%pl%te) == active_nCv .and. size(st%pl%ti) == active_nCv &
            .and. size(st%pl%po) == active_nCv .and. &
            size(st%pl%tn) == active_nCv .and. &
            size(st%pl%kt) == active_nCv .and. size(st%pl%zt) == active_nCv
        if (.not. valid) return
!       Interior cells are 1..nCi, guard cells nCi+1..nCv (b2us_map).
        active_nCi = mpg%nCi
        if (active_nCi .lt. 1 .or. active_nCi .gt. active_nCv) &
            active_nCi = active_nCv

        active_include_po = switch%pot_eq == 1
        active_include_tn = switch%tn_style == 2
        active_include_kt = switch%solve_keps .gt. 0
        active_include_zt = switch%solve_keps .gt. 1
        n_active_species = active_last_species - active_first_species + 1
        n_block = 2 * n_active_species + 2
        if (active_include_po) n_block = n_block + 1
        if (active_include_tn) n_block = n_block + 1
        if (active_include_kt) n_block = n_block + 1
        if (active_include_zt) n_block = n_block + 1
        active_state_size = n_block * active_nCv
        valid = active_state_size .gt. 0
        if (.not. valid) return

        if (allocated(interior_mask)) then
            if (size(interior_mask) .ne. active_state_size) &
                deallocate(interior_mask)
        end if
        if (.not. allocated(interior_mask)) then
            allocate(interior_mask(active_state_size), stat=allocation_status)
            if (allocation_status .ne. 0) then
                valid = .false.
                return
            end if
        end if
        interior_mask = .false.
        do k = 1, n_block
            first = (k - 1) * active_nCv + 1
            last = first + active_nCi - 1
            interior_mask(first:last) = .true.
        end do
    end subroutine configure_active_state
#endif

#if defined(USE_SUNDIALS)
!   The same projection of na onto [na_min, na_max] that b2news_m applies on
!   entry to its density solve (same species range, same
!   use_min_na_numerics branch), applied before the scales and res_first are
!   taken so that both describe a state the map can actually return.
    subroutine clamp_state_na(switch, mpg, st)
        type(switches), intent(in) :: switch
        type(mapping), intent(in) :: mpg
        type(B2State), intent(inout) :: st
        integer :: is, iCv
        real(R8) :: na_min

        do is = active_first_species, active_last_species
            do iCv = 1, active_nCv
                if (switch%use_min_na_numerics .eq. 0) then
                    na_min = switch%b2mndr_na_min
                else
                    na_min = min_na(is, mpg%cvReg(iCv))
                end if
                st%pl%na(iCv, is) = min(max(st%pl%na(iCv, is), na_min), &
                    switch%b2mndr_na_max)
            end do
        end do
    end subroutine clamp_state_na
#endif

#if defined(USE_SUNDIALS)
!   Per-field scales from the interior cells of the entry state, so that a
!   guard cell holding a boundary value far outside the interior range
!   cannot set the units of every unknown of that field.
    subroutine initialize_scale(st, valid)
        type(B2State), intent(in) :: st
        logical, intent(out) :: valid
        integer :: n_active_species, is, scale_index, allocation_status

        valid = state_plasma_is_valid(st)
        if (.not. valid) return

        n_active_species = active_last_species - active_first_species + 1
        if (allocated(scale_na)) then
            if (size(scale_na) .ne. n_active_species) deallocate(scale_na)
        end if
        if (allocated(scale_ua)) then
            if (size(scale_ua) .ne. n_active_species) deallocate(scale_ua)
        end if
        if (.not. allocated(scale_na)) then
            allocate(scale_na(n_active_species), stat=allocation_status)
            if (allocation_status .ne. 0) then
                valid = .false.
                return
            end if
        end if
        if (.not. allocated(scale_ua)) then
            allocate(scale_ua(n_active_species), stat=allocation_status)
            if (allocation_status .ne. 0) then
                valid = .false.
                return
            end if
        end if
        scale_index = 0
        do is = active_first_species, active_last_species
            scale_index = scale_index + 1
            scale_na(scale_index) = get_field_scale(st%pl%na(1:active_nCi, is))
            scale_ua(scale_index) = get_field_scale(st%pl%ua(1:active_nCi, is))
        end do
        scale_po = get_field_scale(st%pl%po(1:active_nCi))
        scale_te = get_field_scale(st%pl%te(1:active_nCi))
        scale_ti = get_field_scale(st%pl%ti(1:active_nCi))
        if (active_include_tn) scale_tn = get_field_scale(st%pl%tn(1:active_nCi))
        if (active_include_kt) scale_kt = get_field_scale(st%pl%kt(1:active_nCi))
        if (active_include_zt) scale_zt = get_field_scale(st%pl%zt(1:active_nCi))
    end subroutine initialize_scale
#endif

#if defined(USE_SUNDIALS)
    subroutine ensure_active_storage(nscxmax, valid)
        integer, intent(in) :: nscxmax
        logical, intent(out) :: valid
        integer :: allocation_status

        valid = .false.
        if (allocated(active_iscx)) then
            if (lbound(active_iscx, 1) .ne. 0 .or. &
                size(active_iscx) .ne. nscxmax) deallocate(active_iscx)
        end if
        if (.not. allocated(active_iscx)) then
            allocate(active_iscx(0:nscxmax-1), stat=allocation_status)
            if (allocation_status .ne. 0) return
        end if
        call ensure_real_buffer(unpack_buffer, active_state_size, valid)
        if (.not. valid) return
        call ensure_c_double_buffer(initial_state_value, active_state_size, &
            valid)
        if (.not. valid) return
        call ensure_c_double_buffer(best_output_value, active_state_size, &
            valid)
        if (.not. valid) return
        call ensure_c_double_buffer(projected_input_value, &
            active_state_size, valid)
    end subroutine ensure_active_storage
#endif

#if defined(USE_SUNDIALS)
    subroutine ensure_workspace(state_size, error_flag, valid)
        integer, intent(in) :: state_size
        integer, intent(out) :: error_flag
        logical, intent(out) :: valid
        integer(c_int) :: flag

        error_flag = 0
        valid = .false.
        if (.not. c_associated(sundials_context)) then
            flag = FSUNContext_Create(SUN_COMM_NULL, sundials_context)
            if (flag .ne. 0_c_int) then
                error_flag = flag
                return
            end if
        end if

        if (workspace_state_size .ne. state_size) then
            call destroy_solver_objects(.false.)
            workspace_state_size = 0
        end if
        if (.not. associated(workspace_state_vector)) then
            workspace_state_vector => FN_VNew_Serial(int(state_size, &
                c_int64_t), sundials_context)
            workspace_u_scale_vector => FN_VNew_Serial(int(state_size, &
                c_int64_t), sundials_context)
            workspace_f_scale_vector => FN_VNew_Serial(int(state_size, &
                c_int64_t), sundials_context)
            if (.not. associated(workspace_state_vector) .or. &
                .not. associated(workspace_u_scale_vector) .or. &
                .not. associated(workspace_f_scale_vector)) then
                error_flag = -1
                call destroy_solver_objects(.false.)
                return
            end if
            workspace_state_size = state_size
        end if
        if (.not. c_associated(kinsol_memory)) then
            kinsol_memory = FKINCreate(sundials_context)
            if (.not. c_associated(kinsol_memory)) then
                error_flag = -1
                return
            end if
        end if
        valid = .true.
    end subroutine ensure_workspace
#endif

#if defined(USE_SUNDIALS)
    subroutine ensure_real_buffer(buffer, required_size, valid)
        real(R8), allocatable, intent(inout) :: buffer(:)
        integer, intent(in) :: required_size
        logical, intent(out) :: valid
        integer :: allocation_status

        if (allocated(buffer)) then
            if (size(buffer) .ne. required_size) deallocate(buffer)
        end if
        if (.not. allocated(buffer)) then
            allocate(buffer(required_size), stat=allocation_status)
            if (allocation_status .ne. 0) then
                valid = .false.
                return
            end if
        end if
        valid = .true.
    end subroutine ensure_real_buffer
#endif

#if defined(USE_SUNDIALS)
    subroutine ensure_c_double_buffer(buffer, required_size, valid)
        real(c_double), allocatable, intent(inout) :: buffer(:)
        integer, intent(in) :: required_size
        logical, intent(out) :: valid
        integer :: allocation_status

        if (allocated(buffer)) then
            if (size(buffer) .ne. required_size) deallocate(buffer)
        end if
        if (.not. allocated(buffer)) then
            allocate(buffer(required_size), stat=allocation_status)
            if (allocation_status .ne. 0) then
                valid = .false.
                return
            end if
        end if
        valid = .true.
    end subroutine ensure_c_double_buffer
#endif

#if defined(USE_SUNDIALS)
    logical function state_is_valid(st) result(valid)
        type(B2State), intent(in) :: st

        valid = all(ieee_is_finite(st%pl%na(:, &
            active_first_species:active_last_species))) .and. &
                all(ieee_is_finite(st%pl%ua(:, &
            active_first_species:active_last_species))) .and. &
                all(ieee_is_finite(st%pl%te)) .and. &
                all(ieee_is_finite(st%pl%ti))
        if (active_include_po) valid = valid .and. &
            all(ieee_is_finite(st%pl%po))
        if (active_include_tn) valid = valid .and. &
            all(ieee_is_finite(st%pl%tn))
        if (active_include_kt) valid = valid .and. &
            all(ieee_is_finite(st%pl%kt))
        if (active_include_zt) valid = valid .and. &
            all(ieee_is_finite(st%pl%zt))
        valid = valid .and. &
            all(st%pl%na(:, active_first_species:active_last_species) .gt. &
                0.0_R8) .and. &
            all(st%pl%te .gt. 0.0_R8) .and. all(st%pl%ti .gt. 0.0_R8) .and. &
            all(st%pl%tn .gt. 0.0_R8)
        if (active_include_kt) valid = valid .and. all(st%pl%kt .ge. 0.0_R8)
        if (active_include_zt) valid = valid .and. all(st%pl%zt .ge. 0.0_R8)
    end function state_is_valid
#endif

#if defined(USE_SUNDIALS)
    logical function state_plasma_is_valid(st) result(valid)
        type(B2State), intent(in) :: st

        valid = state_is_valid(st) .and. &
            all(ieee_is_finite(st%pl%na)) .and. &
            all(ieee_is_finite(st%pl%ua)) .and. &
            all(ieee_is_finite(st%pl%po)) .and. &
            all(ieee_is_finite(st%pl%tn)) .and. &
            all(ieee_is_finite(st%pl%kt)) .and. &
            all(ieee_is_finite(st%pl%zt)) .and. &
            all(st%pl%na .gt. 0.0_R8) .and. all(st%pl%tn .gt. 0.0_R8)
        valid = valid .and. all(st%pl%kt .ge. 0.0_R8) .and. &
            all(st%pl%zt .ge. 0.0_R8)
    end function state_plasma_is_valid
#endif

#if defined(USE_SUNDIALS)
    real(R8) function get_field_scale(buffer) result(scale)
        real(R8), intent(in) :: buffer(:)

        scale = maxval(abs(buffer))
        if (scale .le. 0.0_R8) scale = 1.0_R8
    end function get_field_scale
#endif

#if defined(USE_SUNDIALS)
    subroutine assign_field(field, physical_value, first)
        real(R8), intent(out) :: field(:)
        real(R8), intent(in) :: physical_value(:)
        integer, intent(inout) :: first
        integer :: last

        last = first + size(field) - 1
        field = physical_value(first:last)
        first = last + 1
    end subroutine assign_field
#endif

#if defined(USE_SUNDIALS)
    subroutine unpack_positive_value(buffer, scale, count, physical_value, &
        first, valid)
        real(c_double), intent(in) :: buffer(:)
        real(R8), intent(in) :: scale
        integer, intent(in) :: count
        real(R8), intent(inout) :: physical_value(:)
        integer, intent(inout) :: first
        logical, intent(out) :: valid
        integer :: last
        real(R8) :: lower_limit, upper_limit

        last = first + count - 1
        lower_limit = log(tiny(1.0_R8))
        upper_limit = min(log(huge(1.0_R8)), log(huge(1.0_R8)) - log(scale))
        valid = all(real(buffer(first:last), R8) .ge. lower_limit) .and. &
            all(real(buffer(first:last), R8) .le. upper_limit)
        if (.not. valid) return
        physical_value(first:last) = scale * exp(real(buffer(first:last), R8))
        valid = all(ieee_is_finite(physical_value(first:last))) .and. &
            all(physical_value(first:last) .gt. 0.0_R8)
        if (valid) first = last + 1
    end subroutine unpack_positive_value
#endif

#if defined(USE_SUNDIALS)
    subroutine unpack_signed_value(buffer, scale, count, physical_value, &
        first, valid)
        real(c_double), intent(in) :: buffer(:)
        real(R8), intent(in) :: scale
        integer, intent(in) :: count
        real(R8), intent(inout) :: physical_value(:)
        integer, intent(inout) :: first
        logical, intent(out) :: valid
        integer :: last
        real(R8) :: upper_limit

        last = first + count - 1
        upper_limit = min(log(huge(1.0_R8)), log(huge(1.0_R8)) - log(scale))
        valid = all(abs(real(buffer(first:last), R8)) .le. upper_limit)
        if (.not. valid) return
        physical_value(first:last) = scale * sinh(real(buffer(first:last), R8))
        valid = all(ieee_is_finite(physical_value(first:last)))
        if (valid) first = last + 1
    end subroutine unpack_signed_value
#endif

#if defined(USE_SUNDIALS)
    subroutine unpack_nonnegative_value(buffer, scale, count, &
        physical_value, first, valid)
        real(c_double), intent(in) :: buffer(:)
        real(R8), intent(in) :: scale
        integer, intent(in) :: count
        real(R8), intent(inout) :: physical_value(:)
        integer, intent(inout) :: first
        logical, intent(out) :: valid
        integer :: saved_first

        saved_first = first
        call unpack_signed_value(buffer, scale, count, physical_value, first, &
            valid)
        if (valid) valid = all(physical_value(saved_first:first-1) .ge. 0.0_R8)
        if (.not. valid) first = saved_first
    end subroutine unpack_nonnegative_value
#endif

#if defined(USE_SUNDIALS)
    subroutine unpack_state(buffer, st, valid)
        real(c_double), intent(in) :: buffer(:)
        type(B2State), intent(inout) :: st
        logical, intent(out) :: valid
        integer :: first, last, is, scale_index

        valid = size(buffer) == active_state_size .and. &
            allocated(unpack_buffer) .and. all(ieee_is_finite(buffer))
        if (.not. valid) return

        first = 1
        scale_index = 0
        do is = active_first_species, active_last_species
            scale_index = scale_index + 1
            call unpack_positive_value(buffer, scale_na(scale_index), &
                size(st%pl%na, 1), unpack_buffer, first, valid)
          if (.not. valid) return
        end do
        scale_index = 0
        do is = active_first_species, active_last_species
            scale_index = scale_index + 1
            call unpack_signed_value(buffer, scale_ua(scale_index), &
                size(st%pl%ua, 1), unpack_buffer, first, valid)
            if (.not. valid) return
        end do
        if (active_include_po) then
            call unpack_signed_value(buffer, scale_po, size(st%pl%po), &
                unpack_buffer, first, valid)
            if (.not. valid) return
        end if
        call unpack_positive_value(buffer, scale_te, size(st%pl%te), &
            unpack_buffer, first, valid)
        if (.not. valid) return
        call unpack_positive_value(buffer, scale_ti, size(st%pl%ti), &
            unpack_buffer, first, valid)
        if (.not. valid) return
        if (active_include_tn) then
            call unpack_positive_value(buffer, scale_tn, size(st%pl%tn), &
                unpack_buffer, first, valid)
            if (.not. valid) return
        end if
        if (active_include_kt) then
            call unpack_nonnegative_value(buffer, scale_kt, size(st%pl%kt), &
                unpack_buffer, first, valid)
            if (.not. valid) return
        end if
        if (active_include_zt) then
            call unpack_nonnegative_value(buffer, scale_zt, size(st%pl%zt), &
                unpack_buffer, first, valid)
        end if
        if (.not. valid .or. first .ne. size(buffer) + 1) then
            valid = .false.
            return
        end if

        first = 1
        do is = active_first_species, active_last_species
            last = first + size(st%pl%na, 1) - 1
            st%pl%na(:, is) = unpack_buffer(first:last)
            first = last + 1
        end do
        do is = active_first_species, active_last_species
            last = first + size(st%pl%ua, 1) - 1
            st%pl%ua(:, is) = unpack_buffer(first:last)
            first = last + 1
        end do
        if (active_include_po) &
            call assign_field(st%pl%po, unpack_buffer, first)
        call assign_field(st%pl%te, unpack_buffer, first)
        call assign_field(st%pl%ti, unpack_buffer, first)
        if (active_include_tn) &
            call assign_field(st%pl%tn, unpack_buffer, first)
        if (active_include_kt) &
            call assign_field(st%pl%kt, unpack_buffer, first)
        if (active_include_zt) &
            call assign_field(st%pl%zt, unpack_buffer, first)
        if (.not. active_include_tn) st%pl%tn = st%pl%ti
        valid = state_is_valid(st)
    end subroutine unpack_state
#endif

#if defined(USE_SUNDIALS)
    subroutine pack_positive(field, scale, buffer, first)
        real(R8), intent(in) :: field(:), scale
        real(c_double), intent(inout) :: buffer(:)
        integer, intent(inout) :: first
        integer :: last

        last = first + size(field) - 1
        buffer(first:last) = real(log(field) - log(scale), c_double)
        first = last + 1
    end subroutine pack_positive
#endif

#if defined(USE_SUNDIALS)
    subroutine pack_signed(field, scale, buffer, first)
        real(R8), intent(in) :: field(:), scale
        real(c_double), intent(inout) :: buffer(:)
        integer, intent(inout) :: first
        integer :: last

        last = first + size(field) - 1
        buffer(first:last) = real(asinh(field / scale), c_double)
        first = last + 1
    end subroutine pack_signed
#endif

#if defined(USE_SUNDIALS)
    subroutine pack_state(st, buffer, valid)
        type(B2State), intent(in) :: st
        real(c_double), intent(out) :: buffer(:)
        logical, intent(out) :: valid
        integer :: first, is, scale_index

        valid = size(buffer) == active_state_size
        if (.not. valid) return
        valid = state_is_valid(st)
        if (.not. valid) return

        first = 1
        scale_index = 0
        do is = active_first_species, active_last_species
            scale_index = scale_index + 1
            call pack_positive(st%pl%na(:, is), scale_na(scale_index), &
                buffer, first)
        end do
        scale_index = 0
        do is = active_first_species, active_last_species
            scale_index = scale_index + 1
            call pack_signed(st%pl%ua(:, is), scale_ua(scale_index), &
                buffer, first)
        end do
        if (active_include_po) &
            call pack_signed(st%pl%po, scale_po, buffer, first)
        call pack_positive(st%pl%te, scale_te, buffer, first)
        call pack_positive(st%pl%ti, scale_ti, buffer, first)
        if (active_include_tn) &
            call pack_positive(st%pl%tn, scale_tn, buffer, first)
        if (active_include_kt) &
            call pack_signed(st%pl%kt, scale_kt, buffer, first)
        if (active_include_zt) &
            call pack_signed(st%pl%zt, scale_zt, buffer, first)
        valid = first == size(buffer) + 1 .and. &
            all(ieee_is_finite(buffer))
    end subroutine pack_state
#endif

    subroutine b2_sundials_finalize()
#if defined(USE_SUNDIALS)
    call clear_active_state()
    call destroy_solver_objects(.true.)
    if (allocated(active_iscx)) deallocate(active_iscx)
    if (allocated(scale_na)) deallocate(scale_na)
    if (allocated(scale_ua)) deallocate(scale_ua)
    if (allocated(unpack_buffer)) deallocate(unpack_buffer)
    if (allocated(interior_mask)) deallocate(interior_mask)
    if (allocated(initial_state_value)) deallocate(initial_state_value)
    if (allocated(best_output_value)) deallocate(best_output_value)
    if (allocated(projected_input_value)) deallocate(projected_input_value)
#endif
    end subroutine b2_sundials_finalize

#if defined(USE_SUNDIALS)
    subroutine destroy_solver_objects(flag_release_context)
        logical, intent(in) :: flag_release_context
        integer(c_int) :: flag

        if (c_associated(kinsol_memory)) then
            call FKINFree(kinsol_memory)
            kinsol_memory = c_null_ptr
        end if
        if (associated(workspace_f_scale_vector)) then
            call FN_VDestroy(workspace_f_scale_vector)
            nullify(workspace_f_scale_vector)
        end if
        if (associated(workspace_u_scale_vector)) then
            call FN_VDestroy(workspace_u_scale_vector)
            nullify(workspace_u_scale_vector)
        end if
        if (associated(workspace_state_vector)) then
            call FN_VDestroy(workspace_state_vector)
            nullify(workspace_state_vector)
        end if
        workspace_state_size = 0
        if (flag_release_context .and. c_associated(sundials_context)) then
            flag = FSUNContext_Free(sundials_context)
            sundials_context = c_null_ptr
        end if
    end subroutine destroy_solver_objects
#endif
end module b2mod_sundials
