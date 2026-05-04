program b2test
  implicit none

  print *, 'Unit tests for B2.5'

  call test_sundials_context()
  call test_sundials_vector_initialize()
  call test_sundials_vector_set()
end program b2test

subroutine test_sundials_context()
  use, intrinsic :: iso_c_binding, only: c_int, c_ptr
  use fsundials_core_mod, only: SUN_COMM_NULL, FSUNContext_Create

  implicit none

  type(c_ptr) :: ctx
  integer(c_int) :: ierr

  ierr = FSUNContext_Create(SUN_COMM_NULL, ctx)

  print *, "test.SUNDIALS.context: PASS"
end subroutine test_sundials_context

subroutine test_sundials_vector_initialize()
  use b2us_plasma, only: B2State
  use, intrinsic :: iso_c_binding, only: c_int, c_ptr
  use fsundials_core_mod, only: SUN_COMM_NULL, FSUNContext_Create
  use b2us_plasma_sundials, only: StateVector

  implicit none

  type(c_ptr) :: ctx
  integer(c_int) :: ierr
  type(B2State) :: st
  type(StateVector) :: state

  allocate(st%pl%na(1, 1))
  allocate(st%pl%ua(1, 1))

  ierr = FSUNContext_Create(SUN_COMM_NULL, ctx)
  call state%initialize(ctx, st%pl)

  print *, "test.SUNDIALS.vector.initialize: PASS"
end subroutine test_sundials_vector_initialize

subroutine test_sundials_vector_set()
  use b2us_plasma, only: B2State
  use, intrinsic :: iso_c_binding, only: c_int, c_ptr
  use fsundials_core_mod, only: SUN_COMM_NULL, FSUNContext_Create
  use b2us_plasma_sundials, only: StateVector

  implicit none

  type(c_ptr) :: ctx
  integer(c_int) :: ierr
  type(B2State) :: st
  type(StateVector) :: state

  allocate(st%pl%na(1, 1))
  allocate(st%pl%ua(1, 1))

  ierr = FSUNContext_Create(SUN_COMM_NULL, ctx)
  call state%initialize(ctx, st%pl)
  call state%set(1.0d0)

  print *, "test.SUNDIALS.vector.set: PASS"
end subroutine test_sundials_vector_set
