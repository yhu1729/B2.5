program b2fire
  use b2mod_types, only: R8
  use b2mod_ad, only: nncf
  use b2mod_main, only: b2mn_init, b2mn_fin

  implicit none

  real(kind=R8), dimension(nncf) :: J

  call b2mn_init()
  call b2fire_step(J)
  call b2mn_fin()
end program b2fire

subroutine b2fire_step(J)
  use b2mod_types, only: R8
  use b2mod_ad, only: nncf
  use b2mod_main, only: nout, ns
  use b2us_data, only: mpg, geo, state, state_ext, state_avg, switch
  use b2us_plasma, only: getB2PlasmaSnapshot
  use b2mod_trace, only: b2trcs
  use b2mod_constants, only: qe
  use b2mod_driver, only: ismain, ismain0, itim, dtim, ntim, itim_plas, &
    no_solve, tim

  implicit none

  real(kind=R8), dimension(nncf) :: J
  integer :: ierr, nCv, nFc, nVx
  logical :: ok

  ierr = 0

  nCv = mpg%nCv
  nFc = mpg%nFc
  nVx = mpg%nVx

  call getB2PlasmaSnapshot(state%pl, state%dv, state%psnl)
  call getB2PlasmaSnapshot(state%pl, state%dv, state%psnc)
  call b2trcs()

  if (switch%pot_eq .eq. 0) then
    state%pl%po = 0.0d0
  else if (switch%pot_eq .eq. 2) then
    state%pl%po = 3.1d0 * state%pl%te / qe
  end if

  ok = .false.
  do while(.not. ok)
    call b2mndt (nout, nCv, nFc, nVx, ns, ismain, ismain0, state%rt%nscx, &
      state%rt%nscxmax, state%rt%iscx, itim, dtim, ntim, switch, geo, mpg, &
      state, state_ext, state_avg, ierr)
    ok = ierr .eq. 0

    itim = itim + 1
    itim_plas = itim_plas + 1
    if (no_solve .le. 0) then
      tim = tim + dtim
    endif
  end do
end subroutine b2fire_step
