module b2mod_sundials
    use b2mod_types, only: R8
    use b2mod_switches, only: switches
    use b2mod_ad, only: ncall_b2news_
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
        KIN_FP, KIN_SUCCESS, KIN_INITIAL_GUESS_OK, KIN_STEP_LT_STPTOL, &
        KIN_MAXITER_REACHED, &
        FKINSol, &
        FKINFree, FKINCreate, FKINInit, FKINSetNumMaxIters, FKINSetMAA, &
        FKINSetDamping, FKINSetDampingAA, FKINSetReturnNewest, &
        FKINSetFuncNormTol, FKINSetScaledStepTol, FKINGetNumNonlinSolvIters
#endif

    implicit none
    public

    integer, parameter :: B2_KINSOL_FATAL = -1
    integer, parameter :: B2_KINSOL_SUCCESS = 0
    integer, parameter :: B2_KINSOL_FALLBACK = 1
    integer, parameter :: B2_KINSOL_MAXITER = 2
    integer, parameter :: B2_KINSOL_STEP_TOL = 3

#if defined(USE_SUNDIALS)
    logical, save :: solver_active = .false.
    logical, save :: active_include_po = .true.
    logical, save :: active_include_tn = .true.
    logical, save :: active_include_kt = .false.
    logical, save :: active_include_zt = .false.
    logical, save :: last_output_valid = .false.
    integer, save :: n_active_callback = 0
    integer, save :: n_active_sync = 0
    integer, save :: active_ncall_b2news = 0
    integer, save :: active_first_species = 0
    integer, save :: active_last_species = -1
    integer, save :: active_nCv = 0
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
    real(R8), save :: active_residual_first = 0.0_R8
    real(R8), save :: active_residual_min = 0.0_R8
    real(R8), save :: active_residual_last = 0.0_R8
    real(R8), allocatable, save :: scale_na(:), scale_ua(:)
    real(R8), allocatable, save :: unpack_buffer(:)
    real(R8), save :: scale_po = 1.0_R8
    real(R8), save :: scale_te = 1.0_R8
    real(R8), save :: scale_ti = 1.0_R8
    real(R8), save :: scale_tn = 1.0_R8
    real(R8), save :: scale_kt = 1.0_R8
    real(R8), save :: scale_zt = 1.0_R8
    real(c_double), allocatable, save :: initial_state_value(:)
    real(c_double), allocatable, save :: last_output_value(:)
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
        st_ext, st_avg, max_iteration, maa, fnorm_tolerance, step_tolerance, &
        damping, fnorm_rtolerance, trust_radius, iteration_rxf, &
        iteration_sral, &
        n_iteration, n_evaluation, n_sync_evaluation, &
        residual_first, residual_min, residual_last, &
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
        real(R8), intent(in) :: fnorm_tolerance, step_tolerance, damping
        real(R8), intent(in) :: fnorm_rtolerance, trust_radius, iteration_rxf
        integer, intent(in) :: iteration_sral
        integer, intent(out) :: n_iteration, n_evaluation, n_sync_evaluation
        real(R8), intent(out) :: residual_first, residual_min, residual_last
        integer, intent(out) :: status, kinsol_flag

#if defined(USE_SUNDIALS)
        integer(c_int) :: flag, cleanup_flag
        integer(c_long) :: n_nonlinear_iteration(1), maa_value
        real(c_double), pointer :: state_value(:)
        real(R8) :: effective_fnorm_tolerance, effective_step_tolerance
        real(R8) :: saved_rxf
        integer :: saved_sral
        logical :: returned_state_is_current, valid

        n_iteration = 0
        n_evaluation = 0
        n_sync_evaluation = 0
        residual_first = 0.0_R8
        residual_min = 0.0_R8
        residual_last = 0.0_R8
        active_residual_first = 0.0_R8
        active_residual_min = 0.0_R8
        active_residual_last = 0.0_R8
        active_trust_radius = trust_radius
        active_rxf = iteration_rxf
        active_sral = iteration_sral
        active_fnorm_rtol = fnorm_rtolerance
        status = B2_KINSOL_FALLBACK
        kinsol_flag = 0

        if (solver_active) then
            status = B2_KINSOL_FATAL
            kinsol_flag = -1
            return
        end if
        if (max_iteration .lt. 1) return

        call configure_active_state(switch, st, valid)
        if (.not. valid) goto 900
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
        flag = FKINSetDampingAA(kinsol_memory, &
            real(damping, c_double))
        if (flag .ne. KIN_SUCCESS) then
            kinsol_flag = flag
            goto 900
        end if
        flag = FKINSetReturnNewest(kinsol_memory, 1_c_int)
        if (flag .ne. KIN_SUCCESS) then
            kinsol_flag = flag
            goto 900
        end if
        effective_fnorm_tolerance = fnorm_tolerance
        if (effective_fnorm_tolerance .le. 0.0_R8) then
            effective_fnorm_tolerance = epsilon(1.0_R8)**(1.0_R8/3.0_R8)
        end if
        active_fnorm_floor = effective_fnorm_tolerance
        flag = FKINSetFuncNormTol(kinsol_memory, &
            real(effective_fnorm_tolerance, c_double))
        if (flag .ne. KIN_SUCCESS) then
            kinsol_flag = flag
            goto 900
        end if
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

        active_nCv = nCv
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
        active_ncall_b2news = ncall_b2news_
        last_output_valid = .false.
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

        select case (flag)
        case (KIN_SUCCESS, KIN_INITIAL_GUESS_OK, KIN_STEP_LT_STPTOL, &
            KIN_MAXITER_REACHED)
            returned_state_is_current = n_active_callback == 0
            if (n_active_callback .gt. 0 .and. last_output_valid) then
                returned_state_is_current = all(state_value == last_output_value)
            end if
            if (.not. returned_state_is_current) then
                call unpack_state(state_value, st, valid)
                if (.not. valid) then
                    status = B2_KINSOL_FALLBACK
                    goto 900
                end if

                saved_rxf = switch%b2mndt_rxf
                saved_sral = switch%no_b2sral_call
                switch%b2mndt_rxf = 0.0_R8
                if (iteration_sral .eq. 1) switch%no_b2sral_call = 0
                ncall_b2news_ = active_ncall_b2news
                call b2news_m(nCv, nFc, nVx, ns, nscx, iscx, nscxmax, ismain, &
                    ismain0, dtim, switch, geo, mpg, st, st_ext, st_avg, &
                    .false.)
                switch%b2mndt_rxf = saved_rxf
                switch%no_b2sral_call = saved_sral
                n_active_sync = 1
            end if

            select case (flag)
            case (KIN_MAXITER_REACHED)
                status = B2_KINSOL_MAXITER
            case (KIN_STEP_LT_STPTOL)
                status = B2_KINSOL_STEP_TOL
            case default
                status = B2_KINSOL_SUCCESS
            end select
        case default
            status = B2_KINSOL_FALLBACK
        end select

900 continue
        if (solver_active .and. &
            (status == B2_KINSOL_FATAL .or. &
                status == B2_KINSOL_FALLBACK)) then
            call unpack_state(initial_state_value, st, valid)
            if (valid .and. n_active_callback .gt. 0) then
                saved_rxf = switch%b2mndt_rxf
                saved_sral = switch%no_b2sral_call
                switch%b2mndt_rxf = 0.0_R8
                if (iteration_sral .eq. 1) switch%no_b2sral_call = 0
                ncall_b2news_ = active_ncall_b2news
                call b2news_m(nCv, nFc, nVx, ns, nscx, iscx, nscxmax, &
                    ismain, ismain0, dtim, switch, geo, mpg, st, st_ext, &
                    st_avg, .false.)
                switch%b2mndt_rxf = saved_rxf
                switch%no_b2sral_call = saved_sral
                n_active_sync = n_active_sync + 1
            end if
        end if
        residual_first = active_residual_first
        residual_min = active_residual_min
        residual_last = active_residual_last
        if (n_active_callback .gt. 0) then
            ncall_b2news_ = active_ncall_b2news + n_active_callback + &
                n_active_sync
        end if
        n_evaluation = n_active_callback
        n_sync_evaluation = n_active_sync
        call clear_active_state()
#else
        n_iteration = 0
        n_evaluation = 0
        n_sync_evaluation = 0
        residual_first = 0.0_R8
        residual_min = 0.0_R8
        residual_last = 0.0_R8
        status = B2_KINSOL_FATAL
        kinsol_flag = 0
#endif
    end subroutine b2_sundials_solve

#if defined(USE_SUNDIALS)
    integer(c_int) function b2_sundials_rhs_fixed_point(sunvec_in, &
        sunvec_out, user_data) result(flag) bind(C)
        type(N_Vector) :: sunvec_in, sunvec_out
        type(c_ptr), value :: user_data
        real(c_double), pointer :: value_in(:), value_out(:)
        integer(c_int) :: tolerance_flag
        real(R8) :: change, tolerance, saved_rxf
        integer :: saved_sral
        logical :: valid

        flag = -1_c_int
        if (.not. solver_active) return
        value_in => FN_VGetArrayPointer(sunvec_in)
        value_out => FN_VGetArrayPointer(sunvec_out)
        if (.not. associated(value_in) .or. &
            .not. associated(value_out)) return
        if (size(value_in) .ne. active_state_size .or. &
            size(value_out) .ne. active_state_size) return

        if (active_trust_radius .gt. 0.0_R8) then
            if (maxval(abs(real(value_in, R8) - &
                real(initial_state_value, R8))) .gt. &
                active_trust_radius) return
        end if

        call unpack_state(value_in, active_st, valid)
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

        change = maxval(abs(real(value_out, R8) - real(value_in, R8)))
        active_residual_last = change
        if (n_active_callback .eq. 1) then
            active_residual_first = change
            active_residual_min = change
            if (active_fnorm_rtol .gt. 0.0_R8) then
                tolerance = max(active_fnorm_floor, &
                    active_fnorm_rtol * change)
                if (tolerance .gt. 0.0_R8) then
                    tolerance_flag = FKINSetFuncNormTol(kinsol_memory, &
                        real(tolerance, c_double))
                end if
            end if
        else
            active_residual_min = min(active_residual_min, change)
        end if

        last_output_value = value_out
        last_output_valid = .true.
        flag = 0_c_int
    end function b2_sundials_rhs_fixed_point
#endif

#if defined(USE_SUNDIALS)
    subroutine clear_active_state()
        solver_active = .false.
        n_active_callback = 0
        n_active_sync = 0
        last_output_valid = .false.
        nullify(active_switch, active_geo, active_mpg, active_st, &
            active_st_ext, active_st_avg)
    end subroutine clear_active_state
#endif

#if defined(USE_SUNDIALS)
    subroutine configure_active_state(switch, st, valid)
        type(switches), intent(in) :: switch
        type(B2State), intent(in) :: st
        logical, intent(out) :: valid
        integer :: n_active_species

        active_first_species = max(lbound(st%pl%na, 2), switch%nsmin)
        active_last_species = min(ubound(st%pl%na, 2), switch%nsmax - 1)
        valid = active_last_species .ge. active_first_species
        if (.not. valid) return
        valid = lbound(st%pl%ua, 2) .le. active_first_species .and. &
            ubound(st%pl%ua, 2) .ge. active_last_species
        if (.not. valid) return

        active_include_po = switch%pot_eq == 1
        active_include_tn = switch%tn_style == 2
        active_include_kt = switch%solve_keps .gt. 0
        active_include_zt = switch%solve_keps .gt. 1
        n_active_species = active_last_species - active_first_species + 1
        active_state_size = n_active_species * &
            (size(st%pl%na, 1) + size(st%pl%ua, 1)) + &
            size(st%pl%te) + size(st%pl%ti)
        if (active_include_po) active_state_size = active_state_size + &
            size(st%pl%po)
        if (active_include_tn) active_state_size = active_state_size + &
            size(st%pl%tn)
        if (active_include_kt) active_state_size = active_state_size + &
            size(st%pl%kt)
        if (active_include_zt) active_state_size = active_state_size + &
            size(st%pl%zt)
        valid = active_state_size .gt. 0
    end subroutine configure_active_state
#endif

#if defined(USE_SUNDIALS)
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
            scale_na(scale_index) = maxval(abs(st%pl%na(:, is)))
            scale_ua(scale_index) = maxval(abs(st%pl%ua(:, is)))
            if (scale_ua(scale_index) .le. 0.0_R8) scale_ua(scale_index) = 1.0_R8
        end do
        scale_po = get_field_scale(st%pl%po)
        scale_te = get_field_scale(st%pl%te)
        scale_ti = get_field_scale(st%pl%ti)
        if (active_include_tn) scale_tn = get_field_scale(st%pl%tn)
        if (active_include_kt) scale_kt = get_field_scale(st%pl%kt)
        if (active_include_zt) scale_zt = get_field_scale(st%pl%zt)
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
        call ensure_c_double_buffer(last_output_value, active_state_size, &
            valid)
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
        integer :: first, last, is, scale_index

        valid = size(buffer) == active_state_size
        if (.not. valid) return
        valid = state_is_valid(st)
        if (.not. valid) return

        first = 1
        scale_index = 0
        do is = active_first_species, active_last_species
            scale_index = scale_index + 1
            last = first + size(st%pl%na, 1) - 1
            buffer(first:last) = real(log(st%pl%na(:, is)) - &
                log(scale_na(scale_index)), c_double)
            first = last + 1
        end do
        scale_index = 0
        do is = active_first_species, active_last_species
            scale_index = scale_index + 1
            last = first + size(st%pl%ua, 1) - 1
            buffer(first:last) = real(asinh(st%pl%ua(:, is) / &
                scale_ua(scale_index)), c_double)
            first = last + 1
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
    if (allocated(initial_state_value)) deallocate(initial_state_value)
    if (allocated(last_output_value)) deallocate(last_output_value)
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
