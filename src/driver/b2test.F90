program b2test
#ifdef USE_MPI
  use MPI, only: MPI_Init, MPI_Finalize
#endif

  implicit none

  integer :: ierr

  print *, 'Unit tests for B2.5'

#ifdef USE_MPI
  call MPI_Init(ierr)
#endif

#ifdef USE_MPI
  call test_genex_context()
#endif

#ifdef USE_MPI
  call MPI_Finalize(ierr)
#endif
end program b2test

#ifdef USE_MPI
subroutine test_genex_context()
  use MPI, only: MPI_COMM_WORLD
  use dcomm_handler_m, only: dcomm_handler_t

  implicit none

  type(dcomm_handler_t) :: ctx

  write(*, fmt="(a)", advance="no") "test.GENE-X.context:"

  call ctx%initialize(MPI_COMM_WORLD, 1, 1, 1, 1)

  print *, "PASS"
end subroutine test_genex_context
#endif
