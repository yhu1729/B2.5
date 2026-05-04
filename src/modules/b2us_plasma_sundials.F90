module b2us_plasma_sundials
  use, intrinsic :: iso_c_binding, only: &
    c_int, c_int64_t, c_double, c_ptr, c_loc, c_funloc, c_f_pointer
  use fsundials_core_mod, only: &
    N_Vector, N_Vector_Ops, FN_VNewEmpty, FN_VCopyOps, &
    FN_VGetVecAtIndexVectorArray, SUNDIALS_NVEC_CUSTOM
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
    procedure :: set
    procedure :: find_absolute_value
    procedure :: scale
    procedure :: add
    procedure :: multiply
    procedure :: divide
    procedure :: inverse
    procedure :: axpby
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
    op%nvconst = c_funloc(op_set_all)
    op%nvprod = c_funloc(op_multiply)
    op%nvdiv = c_funloc(op_divide)
    op%nvscale = c_funloc(op_scale)
    op%nvabs = c_funloc(op_absolute_value)
    op%nvinv = c_funloc(op_inverse)
    op%nvdotprod = c_funloc(op_dot_product)
    op%nvmaxnorm = c_funloc(op_max_norm)
    op%nvmin = c_funloc(op_min)
    op%nvwl2norm = c_funloc(op_weighted_l2_norm)
    op%nvl1norm = c_funloc(op_l1_norm)
    op%nvconstrmask = c_funloc(op_constraint_mask)
    op%nvminquotient = c_funloc(op_min_quotient)
    op%nvlinearcombination = c_funloc(op_linear_combination)
    op%nvdotprodmulti = c_funloc(op_dot_product_multiply)
  end subroutine initialize

  ! this = c
  subroutine set(this, c)
    class(StateVector), intent(inout) :: this
    real(kind=8), intent(in) :: c

    this%plasma%na = c
    this%plasma%ua = c
  end subroutine set

  ! this = abs(x)
  subroutine find_absolute_value(this, x)
    class(StateVector), intent(inout) :: this
    class(StateVector), intent(in) :: x

    this%plasma%na = abs(x%plasma%na)
    this%plasma%ua = abs(x%plasma%ua)
  end subroutine find_absolute_value

  ! this = alpha * x
  subroutine scale(this, alpha, x)
    class(StateVector), intent(inout) :: this
    real(kind=8), intent(in) :: alpha
    class(StateVector), intent(in) :: x

    this%plasma%na = alpha * x%plasma%na
    this%plasma%ua = alpha * x%plasma%ua
  end subroutine scale

  ! this = this + alpha * x
  subroutine add(this, alpha, x)
    class(StateVector), intent(inout) :: this
    real(kind=8), intent(in) :: alpha
    class(StateVector), intent(in) :: x

    this%plasma%na = this%plasma%na + alpha * x%plasma%na
    this%plasma%ua = this%plasma%ua + alpha * x%plasma%ua
  end subroutine add

  ! this = x * y
  subroutine multiply(this, x, y)
    class(StateVector), intent(inout) :: this
    class(StateVector), intent(in) :: x
    class(StateVector), intent(in) :: y

    this%plasma%na = x%plasma%na * y%plasma%na
    this%plasma%ua = x%plasma%ua * y%plasma%ua
  end subroutine multiply

  ! this = x / y
  subroutine divide(this, x, y)
    class(StateVector), intent(inout) :: this
    class(StateVector), intent(in) :: x
    class(StateVector), intent(in) :: y

    this%plasma%na = x%plasma%na / y%plasma%na
    this%plasma%ua = x%plasma%ua / y%plasma%ua
  end subroutine divide

  ! this = 1 / x
  subroutine inverse(this, x)
    class(StateVector), intent(inout) :: this
    class(StateVector), intent(in) :: x

    this%plasma%na = 1.0 / x%plasma%na
    this%plasma%ua = 1.0 / x%plasma%ua
  end subroutine inverse

  ! this = s * this + alpha * x + beta * y
  subroutine axpby(this, s, alpha, x, beta, y)
    class(StateVector), intent(inout) :: this
    real(c_double), intent(in) :: s
    real(c_double), intent(in) :: alpha
    class(StateVector), intent(in) :: x
    real(c_double), intent(in) :: beta
    class(StateVector), intent(in) :: y

    this%plasma%na = s * this%plasma%na &
                   + alpha * x%plasma%na &
                   + beta * y%plasma%na
    this%plasma%ua = s * this%plasma%ua &
                   + alpha * x%plasma%ua &
                   + beta * y%plasma%ua
  end subroutine axpby

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

  ! Op: N_VConst
  subroutine op_set_all(c, z) bind(C)
    real(c_double), value, intent(in) :: c
    type(N_Vector), intent(out) :: z
    type(StateVector), pointer :: state_z

    state_z => cast_as_state_vector(z)

    call state_z%set(c)
  end subroutine op_set_all

  ! Op: N_VProd
  subroutine op_multiply(x, y, z) bind(C)
    type(N_Vector), intent(in) :: x
    type(N_Vector), intent(in) :: y
    type(N_Vector), intent(out) :: z
    type(StateVector), pointer :: state_x
    type(StateVector), pointer :: state_y
    type(StateVector), pointer :: state_z

    state_x => cast_as_state_vector(x)
    state_y => cast_as_state_vector(y)
    state_z => cast_as_state_vector(z)

    call state_z%multiply(state_x, state_y)
  end subroutine op_multiply

  ! Op: N_VDiv
  subroutine op_divide(x, y, z) bind(C)
    type(N_Vector), intent(in) :: x
    type(N_Vector), intent(in) :: y
    type(N_Vector), intent(out) :: z
    type(StateVector), pointer :: state_x
    type(StateVector), pointer :: state_y
    type(StateVector), pointer :: state_z

    state_x => cast_as_state_vector(x)
    state_y => cast_as_state_vector(y)
    state_z => cast_as_state_vector(z)

    call state_z%divide(state_x, state_y)
  end subroutine op_divide

  ! Op: N_VScale
  subroutine op_scale(c, x, z) bind(C)
    real(c_double), value, intent(in) :: c
    type(N_Vector), intent(in) :: x
    type(N_Vector), intent(out) :: z
    type(StateVector), pointer :: state_x
    type(StateVector), pointer :: state_z

    state_x => cast_as_state_vector(x)
    state_z => cast_as_state_vector(z)

    call state_z%scale(c, state_x)
  end subroutine op_scale

  ! Op: N_VAbs
  subroutine op_absolute_value(x, z) bind(C)
    type(N_Vector), intent(in) :: x
    type(N_Vector), intent(out) :: z
    type(StateVector), pointer :: state_x
    type(StateVector), pointer :: state_z

    state_x => cast_as_state_vector(x)
    state_z => cast_as_state_vector(z)

    call state_z%find_absolute_value(state_x)
  end subroutine op_absolute_value

  ! Op: N_VInv
  subroutine op_inverse(x, z) bind(C)
    type(N_Vector), intent(in) :: x
    type(N_Vector), intent(out) :: z
    type(StateVector), pointer :: state_x
    type(StateVector), pointer :: state_z

    state_x => cast_as_state_vector(x)
    state_z => cast_as_state_vector(z)

    call state_z%inverse(state_x)
  end subroutine op_inverse

  ! Op: N_VDotProd
  function op_dot_product(x, z) result(d) bind(C)
    real(c_double) :: d
    type(N_Vector), intent(in) :: x
    type(N_Vector), intent(in) :: z
    type(StateVector), pointer :: state_x
    type(StateVector), pointer :: state_z
    real(c_double), dimension(2) :: tmp

    state_x => cast_as_state_vector(x)
    state_z => cast_as_state_vector(z)

    tmp(1) = sum(state_x%plasma%na * state_z%plasma%na)
    tmp(2) = sum(state_x%plasma%ua * state_z%plasma%ua)
    d = sum(tmp)
  end function op_dot_product

  ! Op: N_VMaxNorm
  function op_max_norm(x) result(m) bind(C)
    real(c_double) :: m
    type(N_Vector), intent(in) :: x
    type(StateVector), pointer :: state_x
    real(c_double), dimension(2) :: tmp

    state_x => cast_as_state_vector(x)

    tmp(1) = maxval(abs(state_x%plasma%na))
    tmp(2) = maxval(abs(state_x%plasma%ua))
    m = maxval(tmp)
  end function op_max_norm

  ! Op: N_VMin
  function op_min(x) result(m) bind(C)
    real(c_double) :: m
    type(N_Vector), intent(in) :: x
    type(StateVector), pointer :: state_x
    real(c_double), dimension(2) :: tmp

    state_x => cast_as_state_vector(x)

    tmp(1) = minval(state_x%plasma%na)
    tmp(2) = minval(state_x%plasma%ua)
    m = minval(tmp)
  end function op_min

  ! Op: N_VWL2Norm
  function op_weighted_l2_norm(x, w) result(m) bind(C)
    real(c_double) :: m
    type(N_Vector), intent(in) :: x
    type(N_Vector), intent(in) :: w
    type(StateVector), pointer :: state_x
    type(StateVector), pointer :: state_w

    state_x => cast_as_state_vector(x)

    m = sum(  state_x%plasma%na * state_x%plasma%na &
            * state_w%plasma%na * state_w%plasma%na)
    m = m + sum(  state_x%plasma%ua * state_x%plasma%ua &
                * state_w%plasma%ua * state_w%plasma%ua)
    m = sqrt(m)
  end function op_weighted_l2_norm

  ! Op: N_VL1Norm
  function op_l1_norm(x) result(m) bind(C)
    real(c_double) :: m
    type(N_Vector), intent(in) :: x
    type(StateVector), pointer :: state_x

    state_x => cast_as_state_vector(x)

    m = sum(abs(state_x%plasma%na))
    m = m + sum(abs(state_x%plasma%ua))
  end function op_l1_norm

  ! Op: N_VConstrMask
  function op_constraint_mask(c, x, m) result(t) bind(C)
    integer(c_int) :: t
    type(N_Vector), intent(in) :: c
    type(N_Vector), intent(in) :: x
    type(N_Vector), intent(in) :: m
    type(StateVector), pointer :: state_c
    type(StateVector), pointer :: state_x
    type(StateVector), pointer :: state_m
    integer :: i
    integer :: j
    logical, dimension(2):: test

    state_c => cast_as_state_vector(c)
    state_x => cast_as_state_vector(x)
    state_m => cast_as_state_vector(m)

    t = 1
    do j = 1, size(state_x%plasma%na, 2)
      do i = 1, size(state_x%plasma%na, 1)
        state_m%plasma%na(i, j) = 0.0_8
        state_m%plasma%ua(i, j) = 0.0_8

        if ((state_c%plasma%na(i, j) == 0.0_8) .and. (state_c%plasma%ua(i, j) == 0.0_8)) then
          cycle
        end if
        test(1) = ((abs(state_c%plasma%na(i, j)) > 1.5_8 .and. state_x%plasma%na(i, j) * state_c%plasma%na(i, j) <= 0.0_8) .or. &
                   (abs(state_c%plasma%na(i, j)) > 0.5_8 .and. state_x%plasma%na(i, j) * state_c%plasma%na(i, j) < 0.0_8))
        test(2) = ((abs(state_c%plasma%ua(i, j)) > 1.5_8 .and. state_x%plasma%ua(i, j) * state_c%plasma%ua(i, j) <= 0.0_8) .or. &
                   (abs(state_c%plasma%ua(i, j)) > 0.5_8 .and. state_x%plasma%ua(i, j) * state_c%plasma%ua(i, j) < 0.0_8))
        if (all(test)) then
          t = 0
          state_m%plasma%na(i, j) = 1.0_8
          state_m%plasma%ua(i, j) = 1.0_8
        end if
      end do
    end do
  end function op_constraint_mask

  ! Op: N_VMinQuotient
  function op_min_quotient(num, denom) result(minq) bind(C)
    real(c_double) :: minq
    type(N_Vector), intent(in) :: num
    type(N_Vector), intent(in) :: denom
    type(StateVector), pointer :: state_num
    type(StateVector), pointer :: state_denom
    real(c_double), dimension(2) :: tmp

    tmp(1) = minval(state_num%plasma%na / state_denom%plasma%na)
    tmp(2) = minval(state_num%plasma%na / state_denom%plasma%na)
    minq = minval(tmp)
  end function op_min_quotient

  ! Op: N_VLinearCombination
  function op_linear_combination(nv, c, ptr_x, z) result(ierr) bind(C)
    integer(c_int) :: ierr
    integer(c_int), value, intent(in) :: nv
    real(c_double), dimension(*), intent(in) :: c
    type(c_ptr), value, intent(in) :: ptr_x
    type(N_Vector), intent(out) :: z
    type(StateVector), pointer :: state_z
    integer(c_int) :: k
    type(N_Vector), pointer :: x
    type(StateVector), pointer :: state_x

    state_z => cast_as_state_vector(z)

    x => FN_VGetVecAtIndexVectorArray(ptr_x, 0)
    state_x => cast_as_state_vector(x)
    call state_z%scale(c(1), state_x)

    do k = 2, nv
      x => FN_VGetVecAtIndexVectorArray(ptr_x, k - 1)
      state_x => cast_as_state_vector(x)
      call state_z%add(c(k), state_x)
    end do

    ierr = 0
  end function op_linear_combination

  ! Op: N_VDotProdMulti
  function op_dot_product_multiply(nv, x, ptr_y, d) result(ierr) bind(C)
    integer(c_int) :: ierr
    integer(c_int), value, intent(in) :: nv
    type(N_Vector), intent(in) :: x
    real(c_double), dimension(*), intent(out) :: d
    type(c_ptr), value, intent(in) :: ptr_y
    type(StateVector), pointer :: state_x
    integer(c_int) :: i
    type(N_Vector), pointer :: y
    type(StateVector), pointer :: state_y
    real(c_double), dimension(2) :: tmp

    state_x => cast_as_state_vector(x)
    do i = 1, nv
      y => FN_VGetVecAtIndexVectorArray(ptr_y, i - 1)
      state_y => cast_as_state_vector(y)
      tmp(1) = sum(state_x%plasma%na * state_y%plasma%na)
      tmp(2) = sum(state_x%plasma%ua * state_y%plasma%ua)
      d(i) = sum(tmp)
    end do

    ierr = 0
  end function op_dot_product_multiply
end module b2us_plasma_sundials
