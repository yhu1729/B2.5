program b2test
  implicit none

  write(*, fmt="(a)") "Unit tests for B2.5"

  call test_sundials_context()
  call test_sundials_vector_initialize()
  call test_sundials_vector_set()
end program b2test

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
  use b2us_plasma, only: B2State
  use, intrinsic :: iso_c_binding, only: c_int, c_ptr
  use fsundials_core_mod, only: SUN_COMM_NULL, SUN_SUCCESS, FSUNContext_Create
  use b2us_plasma_sundials, only: StateVector

  implicit none

  type(c_ptr) :: ctx
  integer(c_int) :: ierr
  type(B2State) :: st
  type(StateVector) :: state

  write(*, fmt="(a)", advance="no") "test.SUNDIALS.vector.initialize:"

  allocate(st%pl%na(1, 1))
  allocate(st%pl%ua(1, 1))

  ierr = FSUNContext_Create(SUN_COMM_NULL, ctx)
  if (ierr .ne. SUN_SUCCESS) call abort
  call state%initialize(ctx, st%pl)

  deallocate(st%pl%na)
  deallocate(st%pl%ua)

  print *, "PASS"
end subroutine test_sundials_vector_initialize

subroutine test_sundials_vector_set()
  use b2us_plasma, only: B2State
  use, intrinsic :: iso_c_binding, only: c_int, c_ptr
  use fsundials_core_mod, only: SUN_COMM_NULL, SUN_SUCCESS, FSUNContext_Create
  use b2us_plasma_sundials, only: StateVector

  implicit none

  type(c_ptr) :: ctx
  integer(c_int) :: ierr
  type(B2State) :: st
  type(StateVector) :: state

  write(*, fmt="(a)", advance="no") "test.SUNDIALS.vector.set:"

  allocate(st%pl%na(1, 1))
  allocate(st%pl%ua(1, 1))

  ierr = FSUNContext_Create(SUN_COMM_NULL, ctx)
  if (ierr .ne. SUN_SUCCESS) call abort

  call state%initialize(ctx, st%pl)
  call state%set(1.0d0)

  if (st%pl%na(1, 1) .ne. 1.0d0) call abort
  if (st%pl%ua(1, 1) .ne. 1.0d0) call abort

  deallocate(st%pl%na)
  deallocate(st%pl%ua)

  print *, "PASS"
end subroutine test_sundials_vector_set
