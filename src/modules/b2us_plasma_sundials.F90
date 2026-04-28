module b2us_plasma_sundials
  use, intrinsic :: iso_c_binding, only: &
    c_int, c_int64_t, c_double, c_ptr, c_loc, c_funloc, c_f_pointer
  use fsundials_core_mod, only: &
    N_Vector, N_Vector_Ops, FN_VNewEmpty, FN_VCopyOps, SUNDIALS_NVEC_CUSTOM
  use b2us_plasma, only: &
    B2Plasma

  implicit none

  type StateVector
    type(c_ptr) :: ctx
    type(N_Vector), pointer :: vector
    type(B2Plasma), pointer :: plasma
    integer(c_int64_t) :: length

  contains

    procedure :: initialize
    procedure :: add
    procedure :: axpby
    procedure :: inv
  end type StateVector

contains

  ! Initialize StateVector
  subroutine initialize(this, ctx, plasma)
    class(StateVector), intent(inout), target :: this
    type(c_ptr), intent(in) :: ctx
    type(B2Plasma), target, intent(in) :: plasma
    type(StateVector), pointer :: ptr
    type(N_Vector_Ops), pointer :: op


    this%ctx = ctx
    this%vector => FN_VNewEmpty(this%ctx)
    this%plasma => plasma
    this%length = size(this%plasma%na) + size(this%plasma%ua)

    ptr => this
    this%vector%content = c_loc(ptr)

    call c_f_pointer(this%vector%ops, op)
    op%nvgetvectorid = c_funloc(op_get_vector_id)
    op%nvgetlength = c_funloc(op_get_length)
    op%nvclone = c_funloc(op_clone)
    op%nvdestroy = c_funloc(op_destroy)
    op%nvspace = c_funloc(op_space)
    op%nvlinearsum = c_funloc(op_linear_sum)
    op%nvinv = c_funloc(op_inv)
  end subroutine initialize

  ! this = this + alpha * x
  subroutine add(this, alpha, x)
    class(StateVector), intent(inout) :: this
    real, intent(in) :: alpha
    class(StateVector), intent(in) :: x

    this%plasma%na = this%plasma%na + alpha * x%plasma%na
    this%plasma%ua = this%plasma%ua + alpha * x%plasma%ua
  end subroutine add

  ! this = s * this + alpha * x + beta * y
  subroutine axpby(this, s, alpha, x, beta, y)
    class(StateVector), intent(inout) :: this
    real(c_double), intent(in) :: s
    real(c_double), intent(in) :: alpha
    class(StateVector), intent(in) :: x
    real(c_double), intent(in) :: beta
    class(StateVector), intent(in) :: y

    this%plasma%na = s * this%plasma%na + alpha * x%plasma%na + beta * y%plasma%na
    this%plasma%ua = s * this%plasma%ua + alpha * x%plasma%ua + beta * y%plasma%ua
  end subroutine axpby

  ! z = 1 / this
  subroutine inv(this, z)
    class(StateVector), intent(inout) :: this
    class(StateVector), intent(out) :: z

    z%plasma%na = 1.0 / this%plasma%na
    z%plasma%ua = 1.0 / this%plasma%ua
  end subroutine inv

  ! Cast to StateVector
  module function cast_as_state_vector(x) result(state_x)
    type(StateVector), pointer :: state_x
    type(N_Vector), target, intent(in) :: x

    call c_f_pointer(x%content, state_x)
  end function cast_as_state_vector

  ! Op: N_VGetVectorID
  function op_get_vector_id(w) result(id) bind(C)
    integer(c_int) :: id
    type(N_Vector), intent(in) :: w

    id = SUNDIALS_NVEC_CUSTOM
  end function op_get_vector_id

  ! Op: N_VGetLength
  function op_get_length(v) result(length) bind(C)
    integer(c_int64_t) :: length
    type(N_Vector), intent(in) :: v
    type(StateVector), pointer :: state_v

    state_v => cast_as_state_vector(v)

    length = state_v%length
  end function op_get_length

  ! Op: N_VClone
  function op_clone(w) result(v_ptr) bind(C)
    type(c_ptr) :: v_ptr
    type(N_Vector), intent(inout) :: w
    type(StateVector), pointer :: state_w
    type(N_Vector), pointer :: v
    integer(c_int) :: ierr
    type(StateVector), pointer :: state_v

    state_w => cast_as_state_vector(w)
    v => FN_VNewEmpty(state_w%ctx)
    ierr = FN_VCopyOps(w, v)
    allocate(state_v)
    allocate(state_v%plasma)
    allocate(state_v%plasma%na, mold=state_w%plasma%na)
    allocate(state_v%plasma%ua, mold=state_w%plasma%ua)
    v%content = c_loc(state_v)
    v_ptr = c_loc(v)
  end function op_clone

  ! Op: N_VDestroy
  subroutine op_destroy(v) bind(C)
    type(N_Vector), intent(in) :: v
    type(StateVector), pointer :: state_v

    state_v => cast_as_state_vector(v)

    deallocate(state_v%plasma%na)
    deallocate(state_v%plasma%ua)
    deallocate(state_v%plasma)
  end subroutine op_destroy

  ! Op: N_VSpace
  subroutine op_space(v, length_real, length_integer) bind(C)
    type(N_Vector), intent(in) :: v
    integer(c_int64_t), intent(out) :: length_real(1)
    integer(c_int64_t), intent(out) :: length_integer(1)
    type(StateVector), pointer :: state_v

    state_v => cast_as_state_vector(v)

    length_real(1) = state_v%length
    length_integer(1) = 0
  end subroutine op_space

  ! Op: N_VLinearSum
  subroutine op_linear_sum(a, x, b, y, z) bind(C)
    real(c_double), value, intent(in) :: a
    type(N_Vector), intent(in) :: x
    real(c_double), value, intent(in) :: b
    type(N_Vector), intent(in) :: y
    type(N_Vector), intent(out) :: z
    type(StateVector), pointer :: state_x
    type(StateVector), pointer :: state_y
    type(StateVector), pointer :: state_z

    state_x => cast_as_state_vector(x)
    state_y => cast_as_state_vector(y)
    state_z => cast_as_state_vector(z)

    call state_z%axpby(0.0_8, a, state_x, b, state_y)
  end subroutine op_linear_sum

  ! Op: N_VInv
  subroutine op_inv(x, z) bind(C)
    type(N_Vector), intent(in) :: x
    type(N_Vector), intent(out) :: z
    type(StateVector), pointer :: state_x
    type(StateVector), pointer :: state_z

    state_x => cast_as_state_vector(x)
    state_z => cast_as_state_vector(z)

    call state_x%inv(state_z)
  end subroutine op_inv
end module b2us_plasma_sundials
