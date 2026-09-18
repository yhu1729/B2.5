module b2mod_mms_analytic
!   Time-dependent manufactured solution for temporal verification
!   (b2mndr_use_mms = 3). The exact fields are analytic in space and time;
!   the manufactured sources are built by the discrete-operator method
!     S(t) = d/dt U_exact(t) - L_h(u_exact(t))
!   where L_h is the spatial residual of B2.5 itself, evaluated by the caller
!   (b2mod_sdc: an evaluation-only b2mndt call) and handed to
!   mms_fill_sources. u_exact therefore solves the semi-discrete system
!   exactly and the error measured by mms_write_error is purely temporal
!   (plus the inner nonlinear tolerance).
!
!   Conserved variables match the time terms of b2scdt/b2smdt/b2shdt with
!   tn_style = 0: na, am*mp*na*ua (weight cvVol*cvHz), 1.5*ne*te, 1.5*ni0*ti
!   with ni0 = sum over all species of na (b2xpni). Temperatures are in
!   Joules inside pl; bv_te/bv_ti are in eV (b2stbc_phys multiplies by EV).
    use b2mod_types, only: R8
    use b2mod_constants, only: ev, mp, pi
    use b2mod_b2cmpa, only: am
    use b2us_geo, only: geometry
    use b2us_map, only: mapping
    use b2us_plasma, only: B2Plasma, B2Derivatives
    use b2mod_equation_sources, only: art_sna, art_smo, art_she, art_shi, &
                                      art_shn, art_sch
    use b2mod_boundary_sources, only: bv_na, bv_ua, bv_te, bv_ti, bv_tn, bv_po
    implicit none
    private

    public :: mms_init, mms_exact, mms_exact_dudt, mms_set_boundary_values, &
              mms_fill_sources, mms_set_plasma_state, mms_write_error, &
              mms_initialised

!   parameters (b2mndr_mms_*), read once in mms_init
    real(R8), save :: mms_na0(0:1) = [1.0e16_R8, 1.0e19_R8]  ! background densities [m-3] (neutral, ion)
    real(R8), save :: mms_na_amp = 0.2_R8                    ! relative density amplitude
    real(R8), save :: mms_ua0 = 2.0e4_R8                     ! velocity amplitude [m/s]
    real(R8), save :: mms_te0 = 50.0_R8                      ! electron temperature [eV]
    real(R8), save :: mms_ti0 = 50.0_R8                      ! ion temperature [eV]
    real(R8), save :: mms_t_amp = 0.2_R8                     ! relative temperature amplitude
    real(R8), save :: mms_po_amp = 5.0_R8                    ! potential amplitude [V]
    real(R8), save :: mms_kx = 1.0_R8                        ! poloidal wave number (periods over the box)
    real(R8), save :: mms_ky = 1.0_R8                        ! radial wave number (half periods over the box)
    real(R8), save :: mms_period = 1.0e-3_R8                 ! temporal period [s]
    integer, save :: mms_iout = 0

!   domain box of the interior cells, set in mms_init
    real(R8), save :: box_x0 = 0.0_R8, box_lx = 1.0_R8
    real(R8), save :: box_y0 = 0.0_R8, box_ly = 1.0_R8
    logical, save :: mms_initialised = .false.
    logical, save :: error_file_open = .false.
    integer, parameter :: error_unit = 741
    character(len=*), parameter :: error_filename = 'output/mms_time_error.dat'

    external :: ipgetr, ipgeti, xertst

contains

    subroutine mms_init(nCv, ns, geo, mpg)
        integer, intent(in) :: nCv, ns
        type(geometry), intent(in) :: geo
        type(mapping), intent(in) :: mpg
        real(R8) :: xmin, xmax, ymin, ymax
        integer :: iCv

        if (mms_initialised) return
        call xertst(ns == 2, 'b2mod_mms_analytic: ns must be 2 (D0, D+)')
        call ipgetr('b2mndr_mms_na0_0', mms_na0(0))
        call ipgetr('b2mndr_mms_na0_1', mms_na0(1))
        call ipgetr('b2mndr_mms_na_amp', mms_na_amp)
        call ipgetr('b2mndr_mms_ua0', mms_ua0)
        call ipgetr('b2mndr_mms_te0', mms_te0)
        call ipgetr('b2mndr_mms_ti0', mms_ti0)
        call ipgetr('b2mndr_mms_t_amp', mms_t_amp)
        call ipgetr('b2mndr_mms_po_amp', mms_po_amp)
        call ipgetr('b2mndr_mms_kx', mms_kx)
        call ipgetr('b2mndr_mms_ky', mms_ky)
        call ipgetr('b2mndr_mms_period', mms_period)
        call ipgeti('b2mndr_mms_iout', mms_iout)
        call xertst(all(mms_na0 > 0.0_R8), 'b2mndr_mms_na0_* must be positive')
        call xertst(mms_na_amp >= 0.0_R8 .and. mms_na_amp < 1.0_R8, &
                    'b2mndr_mms_na_amp must be in [0,1)')
        call xertst(mms_te0 > 0.0_R8 .and. mms_ti0 > 0.0_R8, &
                    'b2mndr_mms_te0/ti0 must be positive')
        call xertst(mms_t_amp >= 0.0_R8 .and. mms_t_amp < 1.0_R8, &
                    'b2mndr_mms_t_amp must be in [0,1)')
        call xertst(mms_period > 0.0_R8, 'b2mndr_mms_period must be positive')

!       bounding box of the interior cell centres
        xmin = huge(1.0_R8); xmax = -huge(1.0_R8)
        ymin = huge(1.0_R8); ymax = -huge(1.0_R8)
        do iCv = 1, mpg%nCi
            xmin = min(xmin, geo%cvX(iCv)); xmax = max(xmax, geo%cvX(iCv))
            ymin = min(ymin, geo%cvY(iCv)); ymax = max(ymax, geo%cvY(iCv))
        end do
        box_x0 = xmin; box_lx = max(xmax - xmin, tiny(1.0_R8))
        box_y0 = ymin; box_ly = max(ymax - ymin, tiny(1.0_R8))
        mms_initialised = .true.
        if (nCv < 1) return
        write(*,'(a)') 'b2mod_mms_analytic: analytic MMS (use_mms = 3) active'
        write(*,'(a,2es12.4,a,es12.4)') '  na0 = ', mms_na0, '  amp = ', mms_na_amp
        write(*,'(a,2es12.4,a,es12.4)') '  te0, ti0 [eV] = ', mms_te0, mms_ti0, &
                                        '  amp = ', mms_t_amp
        write(*,'(a,es12.4,a,es12.4)') '  ua0 = ', mms_ua0, '  po_amp = ', mms_po_amp
        write(*,'(a,2f6.2,a,es12.4)') '  kx, ky = ', mms_kx, mms_ky, &
                                      '  period = ', mms_period
        write(*,'(a,4es12.4)') '  box x0, lx, y0, ly = ', box_x0, box_lx, box_y0, box_ly
    end subroutine mms_init

    subroutine shape_functions(x, y, f1, f2)
!       f1: used by na, te, ti;  f2: used by ua, po
        real(R8), intent(in) :: x, y
        real(R8), intent(out) :: f1, f2
        real(R8) :: xi, eta
        xi = (x - box_x0)/box_lx
        eta = (y - box_y0)/box_ly
        f1 = sin(2.0_R8*pi*mms_kx*xi + 0.3_R8)*cos(pi*mms_ky*eta)
        f2 = cos(2.0_R8*pi*mms_kx*xi)*sin(pi*mms_ky*eta + 0.2_R8)
    end subroutine shape_functions

    subroutine time_functions(t, h, dh)
        real(R8), intent(in) :: t
        real(R8), intent(out) :: h, dh
        real(R8) :: w
        w = 2.0_R8*pi/mms_period
        h = sin(w*t)
        dh = w*cos(w*t)
    end subroutine time_functions

    subroutine mms_exact(t, nCv, ns, geo, na, ua, te, ti, po)
!       Exact fields at time t on all cells (guard cells included).
!       te, ti in Joules.
        real(R8), intent(in) :: t
        integer, intent(in) :: nCv, ns
        type(geometry), intent(in) :: geo
        real(R8), intent(out) :: na(nCv,0:ns-1), ua(nCv,0:ns-1), &
                                 te(nCv), ti(nCv), po(nCv)
        real(R8) :: f1, f2, h, dh
        integer :: iCv, is

        call time_functions(t, h, dh)
        do iCv = 1, nCv
            call shape_functions(geo%cvX(iCv), geo%cvY(iCv), f1, f2)
            do is = 0, ns - 1
                na(iCv,is) = mms_na0(is)*(1.0_R8 + mms_na_amp*f1*h)
                ua(iCv,is) = mms_ua0*f2*h
            end do
            te(iCv) = mms_te0*ev*(1.0_R8 + mms_t_amp*f1*h)
            ti(iCv) = mms_ti0*ev*(1.0_R8 + mms_t_amp*f1*h)
            po(iCv) = mms_po_amp*f2*h
        end do
    end subroutine mms_exact

    subroutine mms_exact_dudt(t, nCv, ns, geo, rza, ne_ext, dna, dmom, dEe, dEi)
!       Time derivatives of the conserved variables at time t:
!       dna = d(na)/dt, dmom = d(am*mp*na*ua)/dt, dEe = d(1.5*ne*te)/dt,
!       dEi = d(1.5*ni0*ti)/dt, with ne = sum rza*na + ne_ext (rza assumed
!       constant in time) and ni0 = sum na.
        real(R8), intent(in) :: t
        integer, intent(in) :: nCv, ns
        type(geometry), intent(in) :: geo
        real(R8), intent(in) :: rza(nCv,0:ns-1), ne_ext(nCv)
        real(R8), intent(out) :: dna(nCv,0:ns-1), dmom(nCv,0:ns-1), &
                                 dEe(nCv), dEi(nCv)
        real(R8) :: f1, f2, h, dh, na, ua, dua, te, dte, ti, dti, ne, dne, &
                    ni0, dni0
        integer :: iCv, is

        call time_functions(t, h, dh)
        do iCv = 1, nCv
            call shape_functions(geo%cvX(iCv), geo%cvY(iCv), f1, f2)
            ne = ne_ext(iCv); dne = 0.0_R8
            ni0 = 0.0_R8; dni0 = 0.0_R8
            ua = mms_ua0*f2*h
            dua = mms_ua0*f2*dh
            do is = 0, ns - 1
                na = mms_na0(is)*(1.0_R8 + mms_na_amp*f1*h)
                dna(iCv,is) = mms_na0(is)*mms_na_amp*f1*dh
                dmom(iCv,is) = am(is)*mp*(dna(iCv,is)*ua + na*dua)
                ne = ne + rza(iCv,is)*na
                dne = dne + rza(iCv,is)*dna(iCv,is)
                ni0 = ni0 + na
                dni0 = dni0 + dna(iCv,is)
            end do
            te = mms_te0*ev*(1.0_R8 + mms_t_amp*f1*h)
            dte = mms_te0*ev*mms_t_amp*f1*dh
            ti = mms_ti0*ev*(1.0_R8 + mms_t_amp*f1*h)
            dti = mms_ti0*ev*mms_t_amp*f1*dh
            dEe(iCv) = 1.5_R8*(dne*te + ne*dte)
            dEi(iCv) = 1.5_R8*(dni0*ti + ni0*dti)
        end do
    end subroutine mms_exact_dudt

    subroutine mms_set_boundary_values(t, nCv, ns, geo)
!       Fill the BC type 7 arrays with the exact solution at time t.
        real(R8), intent(in) :: t
        integer, intent(in) :: nCv, ns
        type(geometry), intent(in) :: geo
        real(R8) :: na(nCv,0:ns-1), ua(nCv,0:ns-1), te(nCv), ti(nCv), po(nCv)

        call mms_exact(t, nCv, ns, geo, na, ua, te, ti, po)
        bv_na = na
        bv_ua = ua
        bv_te = te/ev
        bv_ti = ti/ev
        bv_tn = ti/ev
        bv_po = po
    end subroutine mms_set_boundary_values

    subroutine mms_set_plasma_state(t, nCv, ns, geo, rza, ne_ext, pl, dv)
!       pl := u_exact(t) (all cells), derived densities recomputed.
        real(R8), intent(in) :: t
        integer, intent(in) :: nCv, ns
        type(geometry), intent(in) :: geo
        real(R8), intent(in) :: rza(nCv,0:ns-1), ne_ext(nCv)
        type(B2Plasma), intent(inout) :: pl
        type(B2Derivatives), intent(inout) :: dv
        external :: b2xpne, b2xpni, b2xpnn

        call mms_exact(t, nCv, ns, geo, pl%na, pl%ua, pl%te, pl%ti, pl%po)
        pl%tn = pl%ti
        call b2xpni(nCv, ns, pl%na, dv%ni)
        call b2xpnn(nCv, ns, pl%na, dv%nn)
        call b2xpne(nCv, ns, rza, pl%na, ne_ext, dv%ne)
    end subroutine mms_set_plasma_state

    subroutine mms_fill_sources(t, nCv, ns, geo, mpg, rza, ne_ext, &
                                resco0, resmo0, reshe0, reshi0, respo)
!       art_* := d/dt U_exact(t) - L_h(u_exact(t)) over the interior cells,
!       given the volume-integrated spatial residuals res*0 (b2ursc/b2urmo/
!       b2ursd conventions) evaluated at u_exact(t) with art_* = 0.
!       The equations with sign-checked sources (na, te, ti) receive the
!       positive part in component 0 and the negative part as an implicit
!       coefficient evaluated at the exact state, which leaves the residual
!       at u_exact unchanged.
        real(R8), intent(in) :: t
        integer, intent(in) :: nCv, ns
        type(geometry), intent(in) :: geo
        type(mapping), intent(in) :: mpg
        real(R8), intent(in) :: rza(nCv,0:ns-1), ne_ext(nCv)
        real(R8), intent(in) :: resco0(nCv,0:ns-1), resmo0(nCv,0:ns-1), &
                                reshe0(nCv), reshi0(nCv), respo(nCv)
        real(R8) :: na(nCv,0:ns-1), ua(nCv,0:ns-1), te(nCv), ti(nCv), po(nCv)
        real(R8) :: dna(nCv,0:ns-1), dmom(nCv,0:ns-1), dEe(nCv), dEi(nCv)
        real(R8) :: s, ne, ni0
        integer :: iCv, is

        call mms_exact(t, nCv, ns, geo, na, ua, te, ti, po)
        call mms_exact_dudt(t, nCv, ns, geo, rza, ne_ext, dna, dmom, dEe, dEi)
        art_sna = 0.0_R8
        art_smo = 0.0_R8
        art_she = 0.0_R8
        art_shi = 0.0_R8
        art_shn = 0.0_R8
        art_sch = 0.0_R8
        do iCv = 1, mpg%nCi
            ne = ne_ext(iCv)
            ni0 = 0.0_R8
            do is = 0, ns - 1
                ne = ne + rza(iCv,is)*na(iCv,is)
                ni0 = ni0 + na(iCv,is)
!               continuity: sr%sna += cvVol*art_sna
                s = dna(iCv,is) - resco0(iCv,is)/geo%cvVol(iCv)
                art_sna(iCv,0,is) = max(s, 0.0_R8)
                art_sna(iCv,1,is) = min(s, 0.0_R8)/na(iCv,is)
!               momentum: sr%smo += cvVol*cvHz*art_smo and the time term is
!               cvVol*cvHz*d(am*mp*na*ua)/dt (b2smdt), so only the residual
!               carries the weight
                art_smo(iCv,0,is) = dmom(iCv,is) - resmo0(iCv,is)/ &
                                    (geo%cvVol(iCv)*geo%cvHz(iCv))
            end do
!           electron energy: sr%she += cvVol*art_she; implicit slot 3 multiplies ne*te
            s = dEe(iCv) - reshe0(iCv)/geo%cvVol(iCv)
            art_she(iCv,0) = max(s, 0.0_R8)
            art_she(iCv,3) = min(s, 0.0_R8)/(ne*te(iCv))
!           ion energy: slot 3 multiplies ni0*ti
            s = dEi(iCv) - reshi0(iCv)/geo%cvVol(iCv)
            art_shi(iCv,0) = max(s, 0.0_R8)
            art_shi(iCv,3) = min(s, 0.0_R8)/(ni0*ti(iCv))
!           potential (algebraic): sr%sch += cvVol*art_sch
            art_sch(iCv,0) = -respo(iCv)/geo%cvVol(iCv)
        end do
        if (mms_iout /= 0) then
            write(*,'(a,es14.6)') 'mms_fill_sources: t = ', t
            write(*,'(a,2es12.4)') '  max|art_sna| per species = ', &
                (maxval(abs(art_sna(1:mpg%nCi,0,is)) + &
                        abs(art_sna(1:mpg%nCi,1,is)*na(1:mpg%nCi,is))), is=0,ns-1)
            write(*,'(a,2es12.4)') '  max|art_smo| per species = ', &
                (maxval(abs(art_smo(1:mpg%nCi,0,is))), is=0,ns-1)
            write(*,'(a,3es12.4)') '  max|art_she|, |art_shi|, |art_sch| = ', &
                maxval(abs(art_she(1:mpg%nCi,0))+abs(art_she(1:mpg%nCi,3)*ne*te(1:mpg%nCi))), &
                maxval(abs(art_shi(1:mpg%nCi,0))+abs(art_shi(1:mpg%nCi,3)*ni0*ti(1:mpg%nCi))), &
                maxval(abs(art_sch(1:mpg%nCi,0)))
        end if
    end subroutine mms_fill_sources

    subroutine mms_write_error(t, dt, nCv, ns, geo, mpg, pl, &
                               n_solve_sweeps, n_eval_sweeps, sdc_res, n_sdc_sweeps)
!       Append one line per accepted step to output/mms_time_error.dat.
!       Per field: volume-weighted RMS error and RMS error relative to the
!       volume-weighted RMS of the exact field, over interior cells.
        real(R8), intent(in) :: t, dt
        integer, intent(in) :: nCv, ns
        type(geometry), intent(in) :: geo
        type(mapping), intent(in) :: mpg
        type(B2Plasma), intent(in) :: pl
        integer, intent(in) :: n_solve_sweeps, n_eval_sweeps, n_sdc_sweeps
        real(R8), intent(in) :: sdc_res
        real(R8) :: na(nCv,0:ns-1), ua(nCv,0:ns-1), te(nCv), ti(nCv), po(nCv)
        real(R8) :: vals(2, 2*ns + 3)
        integer :: is, k, nCi, ios
        logical :: exists

        nCi = mpg%nCi
        call mms_exact(t, nCv, ns, geo, na, ua, te, ti, po)
        k = 0
        do is = 0, ns - 1
            k = k + 1
            call norms(pl%na(1:nCi,is), na(1:nCi,is), vals(:,k))
        end do
        do is = 0, ns - 1
            k = k + 1
            call norms(pl%ua(1:nCi,is), ua(1:nCi,is), vals(:,k))
        end do
        k = k + 1; call norms(pl%te(1:nCi)/ev, te(1:nCi)/ev, vals(:,k))
        k = k + 1; call norms(pl%ti(1:nCi)/ev, ti(1:nCi)/ev, vals(:,k))
        k = k + 1; call norms(pl%po(1:nCi), po(1:nCi), vals(:,k))

        if (.not. error_file_open) then
            inquire(file=error_filename, exist=exists)
            open(error_unit, file=error_filename, status='unknown', &
                 position='append', action='write', iostat=ios)
            if (ios /= 0) then
!               output/ may not exist yet: fall back to the working directory
                open(error_unit, file='mms_time_error.dat', status='unknown', &
                     position='append', action='write')
                inquire(file='mms_time_error.dat', exist=exists)
            end if
            error_file_open = .true.
            write(error_unit,'(a)') '# b2mod_mms_analytic time-dependent MMS error, one line per step'
            write(error_unit,'(a)') '# columns: t dt rms_na_000 rel_na_000 rms_na_001 rel_na_001'// &
                ' rms_ua_000 rel_ua_000 rms_ua_001 rel_ua_001 rms_te rel_te rms_ti rel_ti'// &
                ' rms_po rel_po sdc_res n_sdc_sweeps n_solve_sweeps n_eval_sweeps'
            write(error_unit,'(a)') '# rms: volume-weighted RMS error over interior cells;'// &
                ' rel: rms / volume-weighted RMS of the exact field; te, ti in eV'
        end if
        write(error_unit,'(2es24.16,14es16.8,es12.4,3i10)') t, dt, &
            (vals(1,k), vals(2,k), k=1,2*ns+3), sdc_res, n_sdc_sweeps, &
            n_solve_sweeps, n_eval_sweeps
        flush(error_unit)

    contains

        subroutine norms(u, uex, out)
            real(R8), intent(in) :: u(:), uex(:)
            real(R8), intent(out) :: out(2)
            real(R8) :: vsum, e2, x2
            vsum = sum(geo%cvVol(1:nCi))
            e2 = sqrt(sum(geo%cvVol(1:nCi)*(u - uex)**2)/vsum)
            x2 = sqrt(sum(geo%cvVol(1:nCi)*uex**2)/vsum)
            out(1) = e2
            if (x2 > 0.0_R8) then
                out(2) = e2/x2
            else
                out(2) = e2
            end if
        end subroutine norms

    end subroutine mms_write_error

end module b2mod_mms_analytic
