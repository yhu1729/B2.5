      subroutine b2news_m (nCv, nFc, nVx, ns, nscx, iscx, nscxmax,
     & ismain, ismain0, dtim, switch, geo, mpg, st, st_ext, st_avg,
     & ramp_call)
         use b2mod_types
         use b2mod_constants
         use b2mod_boundary_namelist
         use b2mod_numerics_namelist
         use b2mod_b2cmpa
         use b2mod_time
         use b2mod_switches
         use b2us_geo
         use b2us_map
         use b2us_plasma
         use b2mod_ma28_for_us, only: restart_ma28_for_us
         use b2mod_facdrift_ExB

#ifndef NAGFOR
         use b2mod_tallies , only:
     &     rrahireg, rranareg, rsahireg, rqradreg, rsanareg, b2bremreg,
     &     rqahereg, rramoreg, b2radreg, rqbrmreg, rsamoreg, b2divue,
     &     b2divua, b2exbe, b2exba, b2visa, b2joule, b2fraa
         use b2mod_transport_namelist, only:
     &     parm_hce, parm_hci, parm_dna, parm_dpa, parm_vla, parm_vsa,
     &     parm_alf, parm_sig
         use b2mod_recycle, only:
     &     int0r, int0l, int3r, int3l, int6r, int6l, int2r, int2l,
     &     int5r, int5l, int1r, int1l, int4r, int4l, fna_mol
         use b2mod_neutrals_namelist, only:
     &     maxw_eff, nstraid, userfluxparm, b2recyc, arcend, targsp,
     &     b2species_start, b2species_end, lstrascl, lsns, msns, gpfc,
     &     nstrai, nstrat, debug_flags, neutrals_count,
     &     neutrals_time_mod, neutrals_time_switch, neut_scl_lim,
     &     natmid, ntrack
         use b2mod_b2cmpt, only:
     &     cfvla, cflim, cfvsa, cfalf, cfdpa, cfsig, cfdna, cfhce, cfhci
         use b2mod_input_profile, only:
     &     ndata, no_pflux, tr_ip_new_files, region_flags, tdata,
     &     tr_ip_elm_count, no_div, addspec, nss, nrr, nkind_coeff,
     &     nkind_data, poloidal_scaling, transport_ip_filename, nscale,
     &     transport_ip_time_switch, transport_ip_time_mod
         use b2us_feedback, only:
     &     fb_target, fb_rescale, fb_current, cum_volrec, fb_const,
     &     charge_frac, saved_fb_actuator, saved_fb_prev, dt_prev,
     &     feedback_namelist_used
         use b2mod_math, only:
     &     cutlo, cutll, b2mod_math_initialised, small_r4_constant
         use b2mod_wall, only: track_species
         use b2mod_zhfrtf, only:
     &     is_i, nal, ia, av_ualpha, gt_ac, gtalc, avm_u,
     &     rho_a_rel, gavm_uc, z_to_m1_ast, w_out, htdp_out, hwdp_out,
     &     rtdpgamma_out, rwdpgamma_out, piu_out
         use b2mod_b2zhco, only:
     &     c_hw_save, c_r_ta, c_r_ta_nofl, c_r_wa, c_r_tb, c_r_tb_nofl,
     &     c_r_wb, c_r_w, tf_ton_mat, amfact
         use b2mod_b2mrel, only:
     &     nucl2s, znnucl_int, s2nucl, amnucl, n2, n3, n4, n5, n6, n8,
     &     n9, n10, n11, n12, n13, n14, n15, n16, mum, rmu, rmd2, mdm,
     &     amnucl_mdf, n2_mdf, n5_mdf, n6_mdf, n8_mdf, n9_mdf, n10_mdf,
     &     n11_mdf, n12_mdf, mum_mdf, rmu_mdf, mdm_mdf
#endif
         use b2mod_ad
         use b2mod_subsys

         implicit none

         intrinsic mod

         integer, intent(in) :: nCv, nFc, nVx, ns
         integer nscx, nscxmax, iscx(0: nscxmax-1), ismain, ismain0
         real(kind=R8) :: dtim
         type(switches), intent(inout) :: switch
         type(geometry), intent(in) :: geo
         type(mapping), intent(inout) :: mpg
         type(B2State), intent(inout) :: st
         type(B2StateExt), intent(inout) :: st_ext
         type(B2Average), intent(in) :: st_avg
         logical, intent(inout) :: ramp_call

         logical :: wrong_flow
         integer :: is, i, ib, iCv, iFc, ireg
         character :: charns*3
         real(kind=R8) :: rxf0, rxf1, rxf2, rxf3, t0
         real(kind=R8) :: zeff(nCv), rz(nCv), wrk0(nCv), zero(nCv),
     &     zerof(nFc, 0:1), wrkf(nFc), fni0(nFc, 0:1)
         logical :: solvedum(0:CVREGMAX), solvedm2(0:CVREGMAX)
         logical :: inverted_ua
         logical :: new_matrix
         real(kind=R8) :: smin, smax, b2sasum
         integer, save :: no_neutr_scl = 0

         external xertst, ipgetr, smin, smax, b2scopy, b2saxpy, b2sasum,
     &     ipgeti
         external b2xvsg, b2xvps, b2xpne, b2xpni, b2xpnn, b2xpfe,
     &     b2xpfi, b2xppb, b2xpnm, b2tral, b2treq, b2tfnb, b2npmo,
     &     b2npco, b2nppo, b2npht, b2nxdp, b2xppz

         ! Start-up
         call subini('b2news_m')
         ! Set internal parameters on first call
         if(ncall_b2news_ .eq. 0) then
            call dev_update_news_solving(b2news_solving, switch)
            call dev_log_news_solving(b2news_solving)
         end if
         call dev_check_parameter(nCv, nFc, ns, ismain, ismain0, geo,
     &     st)

         ! Set rxf[0-3]
         ! The potential equation needs an accurate solution on each
         ! iteration; hence the different setting of rxf2.
         rxf0 = switch%b2mndt_rxf * switch%xfm0
         rxf1 = switch%b2mndt_rxf * switch%xfm1
         rxf2 = 1.0_R8
         if(switch%b2mndt_rxf .eq. 0.0_R8) then
            rxf2 = 0.0_R8
         end if
         rxf3 = switch%b2mndt_rxf * switch%xfm3

         ! Set potential
         zero = 0.0_R8
         zerof = 0.0_R8
         if(switch%pot_eq .eq. 0) then
            st%pl%po = 0.0_R8
         else if(switch%pot_eq .eq. 2) then
            st%pl%po = 3.1_R8 * st%pl%te / qe
         end if

         ! Compute and store factors
         if(ramp_call .or. ncall_b2news_ .eq. 0) then
            call b2news_facdrift_fac_ExB_inc(nFc, st%dv%fac_ExB,
     &        st%dv%facdrift, st%dv%fac_vis, switch)
         end if
         call xertst(0.0_R8 .le. minval(st%dv%facdrift)
     &     .and. 1.0_R8 .ge. maxval(st%dv%facdrift),
     &     'faulty argument facdrift')
         call xertst(0.0_R8 .le. minval(st%dv%fac_ExB)
     &     .and. 1.0_R8 .ge. maxval(st%dv%fac_ExB),
     &     'faulty argument fac_ExB')
         call xertst(0.0_R8 .le. minval(st%dv%fac_vis)
     &     .and. 1.0_R8 .ge. maxval(st%dv%fac_vis),
     &     'faulty argument fac_ExB')

         ! Three passes
         call dev_compute_1(nCv, nFc, nVx, ns, st, st_ext, nscx,
     &     nscxmax, iscx, ismain, switch, geo, mpg, solvedum, solvedm2,
     &     dtim, rxf0, rz, wrk0)
         call dev_compute_2(nCv, nFc, nVx, ns, st, st_ext, st_avg, nscx,
     &     nscxmax, iscx, ismain, ismain0, switch, geo, mpg, dtim, rxf1,
     &     rxf2, rz, wrk0, wrkf, zeff, t0, wrong_flow)
         call dev_compute_3(nCv, nFc, nVx, ns, st, st_ext, ismain,
     &     switch, geo, mpg, solvedum, solvedm2, dtim, fni0,
     &     no_neutr_scl, rxf0, rxf3, wrkf, zero, zerof)

         ncall_b2news_ = ncall_b2news_ + 1
         call subend()
      end subroutine b2news_m

      ! Update solve flags.
      !   flag: Solve flags.
      !   switch: Algorithm switches.
      subroutine dev_update_news_solving(flag, switch)
         use b2mod_switches, only: switches

         implicit none

         logical, intent(out) :: flag(4)
         type(switches), intent(in) :: switch
         integer :: i

         flag(1:4) = switch%no_solve .eq. 0
         do i = 1, 4
            flag(i) = flag(i) .or.
     &      (switch%no_solve .lt. 0 .and.
     &      mod(abs(switch%no_solve), 2**i) .ge. 2**(i - 1))
         end do
      end subroutine dev_update_news_solving

      ! Emit errors when solve flags are wrong.
      !   flag: Solve flags.
      subroutine dev_log_news_solving(flag)
         implicit none

         logical, intent(in) :: flag(4)

         if(.not. flag(1)) then
            write(*, *) 'Not solving parallel momentum equations!'
         end if
         if(.not. flag(2)) then
            write(*, *) 'Not solving continuity equations!'
         end if
         if(.not. flag(3)) then
            write(*,*) 'Not solving potential equation!'
         end if
         if(.not. flag(4)) then
            write(*,*) 'Not solving heat equations!'
         end if
      end subroutine dev_log_news_solving

      ! Check input parameters.
      !   nCv: Number of cells.
      !   nFc: Number of faces.
      !   ns: Number of plasma species
      !   ismain:
      !   ismain0:
      !   geo: Geometry.
      !   st: Plasma state vector.
      subroutine dev_check_parameter(nCv, nFc, ns, ismain, ismain0, geo,
     &   st)
         use b2mod_b2cmpa, only: is_neutral
         use b2mod_ad, only: ncall_b2news_
         use b2us_geo, only: geometry
         use b2mod_types, only: R8
         use b2us_plasma, only: B2State

         implicit none

         integer, intent(in) :: nCv, nFc, ns
         integer, intent(in) :: ismain, ismain0
         type(geometry), intent(in) :: geo
         type(B2State), intent(in) :: st

         real(kind=R8) :: smin, smax

         ! Check scalar parameters
         call xertst(0 .lt. nCv .and. 0 .lt. nFc,
     &   'faulty argument nCv, nFv')
         call xertst(1 .le. ns, 'faulty argument ns')
         call xertst(0 .le. ismain .and. ismain .lt. ns,
     &    'invalid main plasma species index ismain')
         call xertst(.not. is_neutral(ismain),
     &    'invalid main plasma species ismain; must not be neutral')
         call xertst(0 .le. ismain0 .and. ismain0 .lt. ns,
     &    'invalid main neutral species index ismain0')
         call xertst(is_neutral(ismain0) .or. ismain .eq. ismain0,
     &    'invalid main neutral species ismain0; must be neutral')

         ! Check geometry in the in the first 3 calls
         if(ncall_b2news_ .lt. 3) then
            ! Check bb sign
            call b2xvsg(nCv, geo%cvBb(1,3), 1, 'bb3', '.gt.')
            ! Check cvVol and cvQgam
            call b2xvsg(nCv, geo%cvVol, 1, 'vol', '.gt.')
            call xertst(-1.0_R8 .lt. smin(nCv, geo%cvQgam(1, 1), 1)
     &      .and. smax(nCv, geo%cvQgam(1, 1), 1) .lt. 1.0_R8,
     &      'faulty argument range geo%cvQgam(,1)')
            call xertst(0.0_R8 .lt. smin(nCv, geo%cvQgam(1, 0), 1)
     &      .and. smax(nCv, geo%cvQgam(1, 0), 1) .le. 1.0_R8,
     &      'faulty argument range geo%cvQgam(,0)')
            call xertst(-1.0_R8 .lt. smin(nFc, geo%fcQgam(1, 1), 1)
     &      .and. smax(nFc, geo%fcQgam(1, 1), 1) .lt. 1.0_R8,
     &      'faulty argument range geo%fcQgam(,1)')
            call xertst(0.0_R8 .lt. smin(nFc, geo%fcQgam(1, 0), 1)
     &      .and. smax(nFc, geo%fcQgam(1, 0), 1) .le. 1.0_R8,
     &      'faulty argument range geo%fcQgam(,0)')
            ! Check state vector
            call b2xvps(nCv, nFc, ns, st%pl, st%dv)
            call b2xvsg(nCv * ns, st%dv%kinrgy, 1, 'kinrgy', '.ge.')
         end if
      end subroutine dev_check_parameter

      ! Compute 1
      ! Update: ni, nn, ne, ne2
      subroutine dev_compute_1(nCv, nFc, nVx, ns, st, st_ext, nscx,
     &   nscxmax, iscx, ismain, switch, geo, mpg, solvedum, solvedm2,
     &   dtim, rxf0, rz, wrk0)
         use b2us_plasma, only: B2State, B2StateExt
         use b2mod_switches, only: switches
         use b2us_geo, only: geometry
         use b2us_map, only: mapping
         use b2mod_ad, only: CVREGMAX, ncall_b2news_, b2news_solving
         use b2mod_types, only: R8
         use b2mod_numerics_namelist, only: solveco, solvemo,
     &     last_solve_9
         use b2mod_b2cmpa, only: am
         use b2mod_constants, only: mp
         use b2mod_ma28_for_us, only: restart_ma28_for_us

         implicit none

         integer, intent(in) :: nCv, nFc, nVx, ns
         type(B2State), intent(inout) :: st
         type(B2StateExt), intent(inout) :: st_ext
         integer, intent(inout) :: nscx, nscxmax, iscx(0:nscxmax-1),
     &     ismain
         type(switches), intent(inout) :: switch
         type(geometry), intent(in) :: geo
         type(mapping), intent(inout) :: mpg
         logical, intent(inout) :: solvedum(0:CVREGMAX),
     &     solvedm2(0:CVREGMAX)
         real(kind=R8), intent(in) :: dtim
         real(kind=R8), intent(inout) :: rxf0, rz(nCv), wrk0(nCv)

         integer :: is, iCv, ireg
         logical :: new_matrix

         call b2xpni(nCv, ns, st%pl%na, st%dv%ni)
         call b2xpnn(nCv, ns, st%pl%na, st%dv%nn)
         call b2xpne(nCv, ns, st%rt%rza, st%pl%na, st_ext%ne, st%dv%ne)
         call b2xpne(nCv, ns, st%rt%rz2, st%pl%na, st_ext%ne2,
     &     st%dv%ne2)
         call b2tral(nCv, nFc, nVx, ns, nscx, nscxmax, iscx, ismain,
     &     switch, geo, mpg, st%pl, st%dv, st%rt, st_ext, st%co,
     &     st%co_ns)
         call b2tanml(nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &     st%co%csig_an, st%pl%po, st%dv%ne, st%pl%na, st%rt%rza,
     &     st%dv%fchanml_a, st%dv%fchanml)
         call b2tvspa(nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &     st%pl%na, st%pl%ua, st%pl%ti, st%rt%rza, st%co%vsaf_cl,
     &     st%co%vsaf_drho, st%co_ns%vsaf_uAdp_albe,
     &     st%co_ns%vsaf_uBdp_al, st%dv%fac_vis, st%dv%fchvispar_a,
     &     st%dv%fchvispar)
         call b2tiner(nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &     st%pl%na, st%pl%ua, st%dv%facdrift, st%dv%fchinert_a,
     &     st%dv%fchinert)
         call b2tvspr(nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &     st%pl%po, st%pl%ti, st%pl%na, st%rt%rza, st%co%vsa0,
     &     st%dv%fac_vis, st%dv%fchvisper_a, st%dv%fchvisper)
         call b2tvsq(nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &     st%dv%fac_vis, st%dv%ne2, st%rt%rz2, st%dv%lnlam, st%pl%te,
     &     st%pl%ti, st%pl%na, st%dv%ni(:, 1), st%rt%rza,
     &     st%co_ns%vsaf_hAdp_albe, st%co_ns%vsaf_hBdp_al,
     &     st%co%hcix_c, st%co%hci_a, st%co_ns%hci_al_ast,
     &     st%dv%fchvisq_a, st%dv%fchvisq)
         call b2tvskt(nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &     st%pl%po, st%pl%na, st%pl%kt, st%dv%vaecrb, st%co%dna_ExB,
     &     st%dv%fchviskt_a, st%dv%fchviskt)
         if(switch%ion_vlct_restrict .ne. 0 ) then
            call b2xprz(nCv, ns, mp, am, st%pl%na, rz, st_ext)
            call b2xppz(nCv, ns, st%dv%ne, st%pl%na, st%pl%te, st%pl%ti,
     &        st%dv%pz, st_ext)
            wrk0 = sqrt(st%dv%pz/rz) * geo%cvBb(:, 0) / geo%cvBb(:, 3)
            call intface(nCv, nFc, mpg%fcCv, geo%fcVol, wrk0,
     &        st%co%cssb)
         else
            st%co%cssb = 0.0_R8
         end if
         call b2npmo(nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &     ncall_b2news_, b2news_solving(1), 0.0_R8, dtim, st%pl, st%dv,
     &     st%rt, st%co, st%co_ns, st%sr, st%srw, st%psnl, st_ext)
         do is = switch%nsmin, switch%nsmax - 1
            do iCv = 1, nCv
               st%pl%na(iCv, is) = min(max(st%pl%na(iCv, is),
     &           switch%b2mndr_na_min), switch%b2mndr_na_max)
            end do
            call b2xppb(nCv, st%rt%rza(:, is), st%pl%na(:, is),
     &        st%pl%te, st%pl%ti, st%pl%tn, is, st%dv%pa(:, is))
            call b2tfnb(nCv, nFc, nVx, is, ismain, switch, geo, mpg,
     &        st%pl, st%dv, st%co, st%rt, .true.)
            solvedum(0:mpg%nnreg(0)) = solveco(is, 0:mpg%nnreg(0))
            solvedm2(0:mpg%nnreg(0)) = solvemo(is, 0:mpg%nnreg(0))

            new_matrix = .false.
            do ireg = 0, mpg%nnreg(0)
               new_matrix = new_matrix
     &           .or. solveco(is,ireg) .neqv. last_solve_9(ireg)
            end do

            ! No need for "if"
            if(new_matrix .or. .true.) then
               call restart_ma28_for_us
            end if

            call b2npco(nCv, nFc, nVx, mpg%nnreg(0), b2news_solving(2),
     &        solvedum, b2news_solving(1), solvedm2, ncall_b2news_, is,
     &        rxf0, dtim, switch, geo, mpg, st%pl, st%dv, st%sr,
     &        st%psnl)
            last_solve_9(0:mpg%nnreg(0)) = solveco(is, 0:mpg%nnreg(0))

            do iCv = 1, nCv
               st%pl%na(iCv, is) = min(max(st%pl%na(iCv, is),
     &           switch%b2mndr_na_min), switch%b2mndr_na_max)
            end do
            call b2tfnb(nCv, nFc, nVx, is, ismain, switch, geo, mpg,
     &        st%pl, st%dv, st%co, st%rt, .false.)
         end do
         call b2tfed (nCv, nFc, nVx, switch, geo, mpg, st%dv%facdrift,
     &     st%dv%fac_ExB, st%dv%ne, st%pl%te, st%pl%po, st%dv%vedia,
     &     st%dv%veecrb)
      end subroutine dev_compute_1

      ! Compute 2
      ! Update: ne, nn, ni, lnlam
      subroutine dev_compute_2(nCv, nFc, nVx, ns, st, st_ext, st_avg,
     &   nscx, nscxmax, iscx, ismain, ismain0, switch, geo, mpg, dtim,
     &   rxf1, rxf2, rz, wrk0, wrkf, zeff, t0, wrong_flow)
         use b2us_plasma, only: B2State, B2StateExt, B2Average
         use b2mod_switches, only: switches
         use b2us_geo, only: geometry
         use b2us_map, only: mapping
         use b2mod_ad, only: b2news_cutlo, ncall_b2news_, b2news_solving
         use b2mod_types, only: R8
         use b2mod_numerics_namelist, only: last_solve_9, solvepo
         use b2mod_b2cmpa, only: am, zamax, zamin
         use b2mod_constants, only: mp, qe
         use b2mod_boundary_namelist, only: nbc, bcchar
         use b2mod_ma28_for_us, only: restart_ma28_for_us

         implicit none

         integer, intent(in) :: nCv, nFc, nVx, ns
         type(B2State), intent(inout) :: st
         type(B2StateExt), intent(inout) :: st_ext
         type(B2Average), intent(in) :: st_avg
         integer, intent(inout) :: nscx, nscxmax, iscx(0:nscxmax-1),
     &     ismain, ismain0
         type(switches), intent(inout) :: switch
         type(geometry), intent(in) :: geo
         type(mapping), intent(inout) :: mpg
         real(kind=R8), intent(in) :: dtim
         real(kind=R8), intent(inout) :: rxf1, rxf2, rz(nCv), wrk0(nCv),
     &     wrkf(nFc), zeff(nCv), t0
         logical, intent(inout) :: wrong_flow

         integer :: is, i, ib, iFc, ireg
         character :: charns * 3
         logical :: inverted_ua
         logical :: new_matrix
         real(kind=R8) :: b2sasum

         call b2xpne (nCv, ns, st%rt%rza, st%pl%na, st_ext%ne, st%dv%ne)
         call nancheck(nCv, st%dv%ne,
     &   'b2news(ne): after call to b2xpne')
         call b2xpni (nCv, ns, st%pl%na, st%dv%ni)
         call b2xpnn (nCv, ns, st%pl%na, st%dv%nn)
         call b2xpne (nCv, ns, st%rt%rz2, st%pl%na,
     &   st_ext%ne2, st%dv%ne2)
         call b2xppz (nCv, ns, st%dv%ne, st%pl%na, st%pl%te, st%pl%ti,
     &                st%dv%pz, st_ext)
         call b2xpnm (nCv, ns, st%rt%rz2, am, st%pl%na, wrk0)
         call b2treq (nCv, switch, geo%cvVol, st%pl%te, st%pl%ti,
     &   st%dv%ne,
     &     st%dv%ni, wrk0, st%dv%lnlam, st_ext, st%co%ceqp)

         if(switch%b2news_iout .ne. 0) then
            do is = 0,ns-1
               write(charns,'(i3.3)') is
               call my_out_us(70,nCv,0,st%pl%na(1,is),
     &         'b2news__na'//charns)
            end do
            call my_out_us(70,nCv,0,st%dv%ni(1,0),'b2news__ni0')
            call my_out_us(70,nCv,0,st%dv%ni(1,1),'b2news__ni1')
            call my_out_us(70,nCv,0,st%dv%ne,'b2news__ne')
            call my_out_us(70,nCv,0,st%dv%pz,'b2news__pz')
            call my_out_us(70,nCv,0,st%co%ceqp,'b2news__ceqp')
         end if
         call b2tanml (nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &     st%co%csig_an, st%pl%po, st%dv%ne, st%pl%na, st%rt%rza,
     &     st%dv%fchanml_a, st%dv%fchanml)
         call b2tvspa (nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &     st%pl%na, st%pl%ua, st%pl%ti, st%rt%rza, st%co%vsaf_cl,
     &     st%co%vsaf_drho, st%co_ns%vsaf_uAdp_albe,
     &     st%co_ns%vsaf_uBdp_al, st%dv%fac_vis, st%dv%fchvispar_a,
     &     st%dv%fchvispar)
         call b2tiner (nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &     st%pl%na, st%pl%ua, st%dv%facdrift,
     &     st%dv%fchinert_a, st%dv%fchinert)

         do is = switch%nsmin, switch%nsmax-1
            call b2xppb (nCv, st%rt%rza(:,is), st%pl%na(:,is), st%pl%te,
     &          st%pl%ti, st%pl%tn, is, st%dv%pa(:,is))
            call b2tfnb (nCv, nFc, nVx, is, ismain,
     &          switch, geo, mpg, st%pl, st%dv, st%co, st%rt, .true.)
         end do
         call b2tfed (nCv, nFc, nVx, switch, geo, mpg, st%dv%facdrift,
     &       st%dv%fac_ExB, st%dv%ne, st%pl%te, st%pl%po, st%dv%vedia,
     &       st%dv%veecrb)

         if(switch%no_b2sral_call .eq. 0) then
            write(*,*) ' re-compute source '
            call b2sral (nCv, nFc, nVx, ns,
     &        nscx, nscxmax, iscx, ismain, ismain0, dtim,
     &        switch, geo, mpg, st, st_ext, st_avg, wrong_flow, .false.)
            if(wrong_flow) then
               inverted_ua = .true.
               write(*,*) 'b2news_m: wrong_flow returned from b2sral'
               do ib = 1, nbc
                  if(bcchar(ib) .eq. 'A'
     &            .or. bcchar(ib) .eq. 'W'
     &            .or. bcchar(ib) .eq. 'E') then
                     do i = 1, mpg%bcFcP(ib,2)
                        iFc = mpg%bcFc(mpg%bcFcP(ib,1)+i-1)
                        inverted_ua = inverted_ua .and.
     .                   st%dv%fna(iFc,0,ismain)
     &                  *mpg%bcFcOr(mpg%bcFcP(ib,1)+i-1)
     .                    .le.0.0_R8
                     end do
                  end if
               end do
               if(inverted_ua) then
                  write(*,*) 'All target fluxes have the wrong sign'
                  write(*,*)
     .             'Please consider using the b2mndr_inverse_ua switch'
                  call xerrab ('Inverted velocity field!')
               end if
            end if
            call b2srsm (nCv, ns, dtim, switch, geo, mpg,
     &        st%pl%na, st%pl%ua, st%pl%te, st%pl%ti,
     &        st%dv%ne, st%dv%ni, st%dv%nn, st%sr, .false.)
            call b2srdt (nCv, ns, dtim, switch, geo, mpg,
     &        st%psnc%na, st%psnc%ua, st%psnc%te, st%psnc%ti,
     &      st%psnc%tn,
     &        st%psnc%ne, st%psnc%ni, st%psnc%nn, st%psnc%kinrgy,
     &        st%psnc%kt, st%psnc%zt,
     &        st%pl%na, st%pl%ua, st%pl%te, st%pl%ti, st%pl%tn,
     &      st%dv%ne,
     &        st%dv%ni, st%dv%nn, st%dv%kinrgy, st%pl%kt, st%pl%zt,
     &      st%sr,
     &        .false.)
            call b2srst (nCv, ns, switch,
     &        st%pl%na, st%pl%ua, st%pl%te, st%pl%ti, st%pl%tn,
     &      st%pl%po,
     &        st%dv%ne, st%dv%ni, st%dv%nn, st%pl%kt, st%pl%zt, st%sr)
         end if

         wrkf = -st%dv%fch_p(:,0)/(qe*(geo%fcPbs + b2news_cutlo))
         call intcell(nFc,nCv,mpg,mpg%intcellP,wrkf,wrk0)
         do is=0,ns-1
            wrk0 = wrk0 + st%rt%rza(:,is)*st%pl%ua(:,is)*st%pl%na(:,is)
         end do
         wrk0 = wrk0 + st_ext%ue*st_ext%ne
         wrk0 = wrk0/st%dv%ne
         if(switch%b2news_iout.ne.0) then
            call b2xzef (nCv, ns, st%rt%rz2, st%pl%na, st%dv%ne, zeff,
     &       st_ext)
            call b2xpic (nCv, ns, qe, st%rt%rza, geo%cvBb,
     &       st%pl%na, st%pl%ua, zeff, st_ext, st%dv%fch_pi_c)
            call b2xpif(nCv, nFc, switch, geo, mpg, st%dv%fch_pi_c,
     &       st%dv%fch_pi_f)
            call my_out_us(70,nCv,0,st%dv%fch_pi_c,'b2news__fch_pi_c')
            call my_out_us(70,nFc,1,st%dv%fch_pi_f,'b2news__fch_pi_f')
            call my_out_us(70,nCv,0,st%dv%ue,'b2news__ve')
            call my_out_us(70,nCv,0,wrk0,'b2news__ue')
         end if
         st%dv%ue=wrk0
         call b2npmo (nCv, nFc, nVx, ns, ismain,
     &     switch, geo, mpg,
     &     ncall_b2news_, b2news_solving(1), rxf1, dtim,
     &     st%pl, st%dv, st%rt, st%co, st%co_ns, st%sr, st%srw,
     &     st%psnl, st_ext)
         if(switch%pot_eq.eq.0) then
            st%pl%po        = 0.0_R8

            call b2tfch_ (nCv, nFc, nVx, ns, ismain, switch,
     &         geo, mpg, st%pl, st%dv, st%co, st%co_ns, st%rt, st_ext)
         else if(switch%pot_eq.eq.2) then
            st%pl%po        = 3.1_R8*st%pl%te/qe

            call b2tfch_ (nCv, nFc, nVx, ns, ismain, switch,
     &         geo, mpg, st%pl, st%dv, st%co, st%co_ns, st%rt, st_ext)
         else
            if(b2news_solving(3)) then
               call b2nxdp (nCv, mpg, solvepo, qe, st%pl%te, st%psnl%ne,
     &           st%dv%ne, wrk0)
               call b2saxpy (nCv, 1.0_R8, wrk0, 1, st%pl%po, 1)
               if(switch%b2news_iout .ne. 0) then
                  call my_out_us(70,nCv,0,st%pl%po,'b2nxdp_po')
               end if
            end if
            call b2tfch_ (nCv, nFc, nVx, ns, ismain, switch,
     &         geo, mpg, st%pl, st%dv, st%co, st%co_ns, st%rt, st_ext)
            if(b2news_solving(3)) then
               if(switch%fhe_vis_per.eq.0.0_R8.or..true.) then
                  new_matrix = .false.
                  do ireg = 0, mpg%nnreg(0)
                     new_matrix = new_matrix.or.
     &                solvepo(ireg).neqv.last_solve_9(ireg)
                  end do
                  if(new_matrix.or..true.) call restart_ma28_for_us
                  call b2nppo (nCv, nFc, nVx, mpg%nnreg(0),
     &              b2news_solving(3), solvepo, ncall_b2news_, rxf2,
     &            switch,
     &              geo, mpg, st%pl, st%dv, st%sr)
                  last_solve_9(0:mpg%nnreg(0))=solvepo(0:mpg%nnreg(0))
               else
                  call xerrab ('visper not yet available for WG')
#ifdef WG_TODO
                  call b2npp7 (nx, ny, mpg%nnreg(0), b2news_solving(3),
     &             solvepo, ncall_b2news_, rxf2, ne, te, conc, avisper,
     &             bvisper, cvisper, sch, po, fch, respo, corpo)
#endif
               end if
               t0 = b2sasum(nCv,st%dv%respo,1)
               write(*,*) 'respo before solving',t0
               call b2tfch_ (nCv, nFc, nVx, ns, ismain, switch,
     &          geo, mpg, st%pl, st%dv, st%co, st%co_ns, st%rt, st_ext)
               call b2ursd (nCv, nFc, mpg, st%dv%ne, st%pl%po,
     &          st%sr%sch, st%dv%fch, wrk0)
               t0 = b2sasum(nCv,wrk0,1)
               write(*,*) 'respo after solving',t0
            else
               if(switch%get_residuals.eq.1)
     &          call b2ursd (nCv, nFc, mpg,
     &                       st%dv%ne, st%pl%po, st%sr%sch, st%dv%fch,
     &                       st%dv%respo)
            end if
         end if

         if(switch%b2news_iout .ne. 0) then
            call my_out_us(70,nCv,0,st%dv%respo,'b2news__respo')
            call my_out_us(70,nCv,0,st%pl%po,'b2news__po')
         end if

         if(switch%ion_vlct_restrict .ne. 0 ) then
            call b2xprz (nCv, ns, mp, am, st%pl%na, rz, st_ext)
            call b2xppz (nCv, ns, st%dv%ne, st%pl%na, st%pl%te,
     &      st%pl%ti,
     &                   st%dv%pz, st_ext)
            wrk0 = sqrt(st%dv%pz/rz)*geo%cvBb(:,0)/geo%cvBb(:,3)
            call intface(nCv,nFc,mpg%fcCv,geo%fcVol,wrk0,st%co%cssb)
         else
            st%co%cssb = 0.0_R8
         end if
         do is = 0, ns-1
            call b2xppb (nCv, st%rt%rza(:,is), st%pl%na(:,is), st%pl%te,
     &        st%pl%ti, st%pl%tn, is, st%dv%pa(:,is))
#ifdef DBG
            write(*,*) 'is, za(is) ',is, (zamin(is)+zamax(is))/2.0_R8
#endif
            call b2tfnb (nCv, nFc, nVx, is, ismain,
     &       switch, geo, mpg, st%pl, st%dv, st%co, st%rt, .false.)
         end do
      end subroutine dev_compute_2

      ! Compute 3
      ! Update: fne_he, fni_he
      subroutine dev_compute_3(nCv, nFc, nVx, ns, st, st_ext, ismain,
     &   switch, geo, mpg, solvedum, solvedm2, dtim, fni0, no_neutr_scl,
     &   rxf0, rxf3, wrkf, zero, zerof)
         use b2us_plasma, only: B2State, B2StateExt
         use b2mod_switches, only: switches
         use b2us_geo, only: geometry
         use b2us_map, only: mapping
         use b2mod_ad, only: CVREGMAX, b2news_cutlo, ncall_b2news_,
     &     b2news_solving
         use b2mod_types, only: R8
         use b2mod_numerics_namelist, only: solveco, solvemo, solveee,
     &     solveei, solveen, solveet, solvepo, solvemo, solvekt,
     &     solvezt, last_solve_9
         use b2mod_constants, only: mp, qe
         use b2mod_ma28_for_us, only: restart_ma28_for_us

         implicit none

         integer, intent(in) :: nCv, nFc, nVx, ns
         type(B2State), intent(inout) :: st
         type(B2StateExt), intent(inout) :: st_ext
         integer, intent(inout) :: ismain
         type(switches), intent(inout) :: switch
         type(geometry), intent(in) :: geo
         type(mapping), intent(inout) :: mpg
         logical, intent(inout) :: solvedum(0:CVREGMAX),
     &     solvedm2(0:CVREGMAX)
         real(kind=R8), intent(in) :: dtim
         real(kind=R8), intent(inout) :: fni0(nFc,0:1), wrkf(nFc),
     &     zero(nCv), zerof(nFc,0:1)
         integer, intent(inout) :: no_neutr_scl
         real(kind=R8), intent(inout) :: rxf0, rxf3

         integer :: is, iCv, ireg
         logical :: new_matrix

         call b2xpfe (nCv, nFc, ns, st_ext%ns, geo, mpg, qe,
     &     st%rt%rza, st_ext%za, st%dv%fna_he, st_ext%fa,
     &     st%dv%fch_p, st%dv%fne_he)
         call b2xpfe (nCv, nFc, ns, st_ext%ns, geo, mpg, qe,
     &     st%rt%rza, st_ext%za, st%dv%fna_53, st_ext%fa,
     &     st%dv%fch_p, st%dv%fne_53)
         call b2xpfe (nCv, nFc, ns, st_ext%ns, geo, mpg, qe,
     &     st%rt%rza, st_ext%za, st%dv%fna_32nd, st_ext%fa,
     &     st%dv%fch_32, st%dv%fne_32)
         call b2xpfe (nCv, nFc, ns, 0, geo, mpg, qe,
     &     st%rt%rza, zero, st%dv%fna_52nd, zerof,
     &     st%dv%fch_52, st%dv%fne_52)
         if(switch%tn_style.eq.0) then

            call b2xpfi (nFc, ns, st%dv%fna_he, st%dv%fni_he)
            call b2xpfi (nFc, ns, st%dv%fna_32nd, st%dv%fni_32)
            call b2xpfi (nFc, ns, st%dv%fna_52nd, st%dv%fni_52)
         else if(switch%tn_style.eq.1) then

            call b2xpfz (nFc, ns, st%dv%fna_he, st%dv%fni_he)
            call b2xpfz (nFc, ns, st%dv%fna_32nd, st%dv%fni_32)
            call b2xpfz (nFc, ns, st%dv%fna_52nd, st%dv%fni_52)
         else


            call b2xpfi (nFc, ns, st%dv%fna_he, fni0)
            call b2xpfn (nFc, ns, st%dv%fna_he, st%dv%fnn)
            st%dv%fni_he = fni0 - st%dv%fnn
            call b2xpfi (nFc, ns, st%dv%fna_32nd, fni0)
            call b2xpfn (nFc, ns, st%dv%fna_32nd, st%dv%fnn_32)
            st%dv%fni_32 = fni0 - st%dv%fnn_32
            call b2xpfi (nFc, ns, st%dv%fna_52nd, fni0)
            call b2xpfn (nFc, ns, st%dv%fna_52nd, st%dv%fnn_52)
            st%dv%fni_52 = fni0 - st%dv%fnn_52
         end if
         call b2tfed (nCv, nFc, nVx, switch, geo, mpg, st%dv%facdrift,
     &     st%dv%fac_ExB, st%dv%ne, st%pl%te, st%pl%po,
     &     st%dv%vedia, st%dv%veecrb)
         call b2tfhe_ (nCv, nFc, nVx, ns, switch, geo, mpg,
     &     st%pl, st%dv, st%co, st%rt, st_ext, .true.)
         call b2tfhi_(nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &     st%pl, st%dv, st%co, .true.)
         st%dv%floi_noc = st%dv%floi
         st%dv%floe_noc = st%dv%floe
         call b2xpen (nCv, nFc, ns, mpg, switch,
     &     st%dv%fna, st%dv%fna_53, st%dv%fne_53,
     &     st%dv%fch, st%dv%fhe, st%dv%fhi, st%dv%fhn, st%rt%rpt,
     &   st%pl%te,
     &     st%pl%ti, st%pl%tn, st%pl%po,
     &     st%dv%fhm, st%dv%fhp, st%dv%fhj, st%dv%fht)

         if(switch%iout_b2wdat.eq.4) then
            call my_out_us(70,nFc,1,st%dv%fht(1,0),'b2news_m_fht_th')
            call my_out_us(70,nFc,1,st%dv%fht(1,1),'b2news_m_fht_r')
         end if

         call b2tfcc (nCv, nFc, nVx, ns, switch, geo, mpg,
     &     st%pl, st%dv, st%co, st%rt)

         wrkf = -st%dv%fch_p(:,0)/(qe*(geo%fcPbs + b2news_cutlo))
         call intcell(nFc,nCv,mpg,mpg%intcellP,wrkf,st%dv%ue)
         do is=0,ns-1
            st%dv%ue = st%dv%ue +
     &               st%rt%rza(:,is)*st%pl%ua(:,is)*st%pl%na(:,is)
         end do
         st%dv%ue = st%dv%ue + st_ext%ue*st_ext%ne
         st%dv%ue = st%dv%ue/st%dv%ne

         call b2npht (nCv, nFc, nVx, ns, switch, geo, mpg,
     &     ncall_b2news_, ismain, b2news_solving, solveee, solveei,
     &     solveen, solveet, solvepo, solvemo, solvekt, solvezt, rxf3,
     &     dtim, st%pl, st%dv, st%co, st%rt, st%sr, st%srw,
     &     st%psnl, st_ext)
         if(switch%pot_eq.eq.0) st%pl%po=0.0_R8
         if(switch%pot_eq.eq.2) st%pl%po=3.1_R8*st%pl%te/qe

         if(switch%do_2nd_b2npco_call.ne.0) then
            do is=0,ns-1
               do iCv = 1, nCv
                  st%pl%na(iCv,is)=min(max(st%pl%na(iCv,is),
     &              switch%b2mndr_na_min),switch%b2mndr_na_max)
               end do
               call b2xppb (nCv, st%rt%rza(:,is), st%pl%na(:,is),
     &         st%pl%te,
     &           st%pl%ti, st%pl%tn, is, st%dv%pa(:,is))
               call b2tfnb (nCv, nFc, nVx, is, ismain,
     &           switch, geo, mpg, st%pl, st%dv, st%co, st%rt, .false.)
               new_matrix = .false.
               solvedum(0:mpg%nnreg(0))=solveco(is,0:mpg%nnreg(0))
               solvedm2(0:mpg%nnreg(0))=solvemo(is,0:mpg%nnreg(0))
               do ireg = 0, mpg%nnreg(0)
                  new_matrix = new_matrix.or.
     &             solveco(is,ireg).neqv.last_solve_9(ireg)
               end do
               if(new_matrix) call restart_ma28_for_us
               call b2npco (nCv, nFc, nVx, mpg%nnreg(0),
     &           b2news_solving(2), solvedum, b2news_solving(1),
     &         solvedm2,
     &           ncall_b2news_, is, rxf0, dtim, switch, geo, mpg,
     &           st%pl, st%dv, st%sr, st%psnl)
               last_solve_9(0:mpg%nnreg(0))=solveco(is,0:mpg%nnreg(0))
            end do
            call b2tfed (nCv, nFc, nVx, switch, geo, mpg,
     &       st%dv%facdrift, st%dv%fac_ExB, st%dv%ne, st%pl%te,
     &      st%pl%po,
     &       st%dv%vedia, st%dv%veecrb)
         end if

         if(switch%re_eval_prtls_fluxes.eq.1.or.
     &   switch%use_eirene.ne.0) then
            do is=0,ns-1
               do iCv = 1, nCv
                  st%pl%na(iCv,is)=min(max(st%pl%na(iCv,is),
     &              switch%b2mndr_na_min),switch%b2mndr_na_max)
               end do
               call b2xppb (nCv, st%rt%rza(:,is), st%pl%na(:,is),
     &         st%pl%te,
     &           st%pl%ti, st%pl%tn, is, st%dv%pa(:,is))
               call b2tfnb (nCv, nFc, nVx, is, ismain,
     &           switch, geo, mpg, st%pl, st%dv, st%co, st%rt, .false.)
            end do
            call b2tfed (nCv, nFc, nVx, switch, geo, mpg,
     &      st%dv%facdrift,
     &        st%dv%fac_ExB, st%dv%ne, st%pl%te, st%pl%po,
     &        st%dv%vedia, st%dv%veecrb)
         end if
         if(no_neutr_scl .eq. 0) then
            write(*,*) 'Calling b2stbr_neutr_scl'
            call b2stbr_neutr_scl(nCv, nFc, ns, mpg, geo, st%pl, st%dv,
     &        switch%use_eirene, switch%ank_mods)
            write(*,*) 'Called b2stbr_neutr_scl'
         end if
      end subroutine dev_compute_3
