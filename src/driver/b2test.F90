program b2test
  implicit none

  write(*, fmt="(a)") "Unit tests for B2.5"

  call test_sundials_context()
  call test_sundials_vector_initialize()
  call test_sundials_vector_set()
  call test_sundials_vector_find_absolute_value()
end program b2test

subroutine setup_fixture(ctx, st, state)
  use, intrinsic :: iso_c_binding, only: c_int, c_ptr
  use b2us_plasma, only: B2State
  use b2us_plasma_sundials, only: StateVector
  use fsundials_core_mod, only: SUN_COMM_NULL, SUN_SUCCESS, FSUNContext_Create

  implicit none

  type(c_ptr), intent(inout) :: ctx
  type(B2State), intent(inout) :: st
  type(StateVector), intent(inout) :: state
  integer(c_int) :: ierr

  allocate(st%pl%na(128, 32))
  allocate(st%pl%ua(128, 32))
  allocate(st%pl%po(128))
  allocate(st%pl%te(128))
  allocate(st%pl%ti(128))
  allocate(st%pl%tn(128))
  allocate(st%pl%kt(128))
  allocate(st%pl%zt(128))

  ierr = FSUNContext_Create(SUN_COMM_NULL, ctx)
  if (ierr .ne. SUN_SUCCESS) call abort

  call state%initialize(ctx, st%pl)
end subroutine setup_fixture

subroutine teardown_fixture(ctx, st)
  use, intrinsic :: iso_c_binding, only: c_ptr
  use b2us_plasma, only: B2State

  implicit none

  type(c_ptr), intent(in) :: ctx
  type(B2State), intent(inout) :: st

  deallocate(st%pl%na)
  deallocate(st%pl%ua)
  deallocate(st%pl%po)
  deallocate(st%pl%te)
  deallocate(st%pl%ti)
  deallocate(st%pl%tn)
  deallocate(st%pl%kt)
  deallocate(st%pl%zt)
end subroutine teardown_fixture

subroutine test_sundials_context()
  use, intrinsic :: iso_c_binding, only: c_int, c_ptr
  use fsundials_core_mod, only: SUN_COMM_NULL, SUN_SUCCESS, FSUNContext_Create

  implicit none

  type(c_ptr) :: ctx
  integer(c_int) :: ierr

  write(*, fmt="(a)", advance="no") "test.SUNDIALS.context:"

  ierr = FSUNContext_Create(SUN_COMM_NULL, ctx)
  if (ierr .ne. SUN_SUCCESS) call abort

  print *, "PASS"
end subroutine test_sundials_context

subroutine test_sundials_vector_initialize()
  use, intrinsic :: iso_c_binding, only: c_ptr
  use b2us_plasma, only: B2State
  use b2us_plasma_sundials, only: StateVector

  implicit none

  type(c_ptr) :: ctx
  type(B2State) :: st
  type(StateVector) :: state

  write(*, fmt="(a)", advance="no") "test.SUNDIALS.vector.initialize:"

  call setup_fixture(ctx, st, state)

  if (state%length .ne. 8960) call abort

  call teardown_fixture(ctx, st)

  print *, "PASS"
end subroutine test_sundials_vector_initialize

subroutine test_sundials_vector_set()
  use, intrinsic :: iso_c_binding, only: c_ptr
  use b2us_plasma, only: B2State
  use b2us_plasma_sundials, only: StateVector

  implicit none

  type(c_ptr) :: ctx
  type(B2State) :: st
  type(StateVector) :: state

  write(*, fmt="(a)", advance="no") "test.SUNDIALS.vector.set:"

  call setup_fixture(ctx, st, state)

  call state%set(1.0d0)

  if (st%pl%na(128, 32) .ne. 1.0d0) call abort
  if (st%pl%ua(128, 32) .ne. 1.0d0) call abort
  if (st%pl%po(128) .ne. 1.0d0) call abort
  if (st%pl%te(128) .ne. 1.0d0) call abort
  if (st%pl%ti(128) .ne. 1.0d0) call abort
  if (st%pl%tn(128) .ne. 1.0d0) call abort
  if (st%pl%kt(128) .ne. 1.0d0) call abort
  if (st%pl%zt(128) .ne. 1.0d0) call abort

  call teardown_fixture(ctx, st)

  print *, "PASS"
end subroutine test_sundials_vector_set

subroutine test_sundials_vector_find_absolute_value()
  use, intrinsic :: iso_c_binding, only: c_ptr
  use b2us_plasma, only: B2State
  use b2us_plasma_sundials, only: StateVector

  implicit none

  type(c_ptr) :: ctx
  type(B2State) :: st
  type(StateVector) :: state
  type(StateVector) :: output

  write(*, fmt="(a)", advance="no") "test.SUNDIALS.vector.find_absolute_value:"

  call setup_fixture(ctx, st, state)

  call output%clone(state)
  call state%set(-1.0d0)
  call output%find_absolute_value(state)
  if (output%plasma%na(128, 32) .ne. 1.0d0) call abort
  if (output%plasma%ua(128, 32) .ne. 1.0d0) call abort
  if (output%plasma%po(128) .ne. 1.0d0) call abort
  if (output%plasma%te(128) .ne. 1.0d0) call abort
  if (output%plasma%ti(128) .ne. 1.0d0) call abort
  if (output%plasma%tn(128) .ne. 1.0d0) call abort
  if (output%plasma%kt(128) .ne. 1.0d0) call abort
  if (output%plasma%zt(128) .ne. 1.0d0) call abort

  call output%destroy()
  call teardown_fixture(ctx, st)

  print *, "PASS"
end subroutine test_sundials_vector_find_absolute_value
