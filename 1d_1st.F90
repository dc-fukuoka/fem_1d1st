! the simplest case, 1d and 1st order(linear) element
! solve d^2u/dx^2 = lambda (0 <= x <= L)
! boundary condition: u(0)     = 0 (Dirichlet boundary condition)
!                     du/dx(L) = 0 (Neumann boundary condition)
! analytic solution: u(x) = lambda/2*x(x-2L)

module vars
  implicit none
  integer,parameter  :: ui = 7
  integer,parameter  :: dp = kind(1.0d0)
  real(dp),parameter :: length = 3.0d0 ! = L
  integer,parameter  :: iter_max = 10000
  real(dp),parameter :: tol = 1.0d-10
  real(dp),parameter :: lambda = 1.0d0
  type elem
     integer :: nelements
     integer :: nnodes
     real(dp),allocatable,dimension(:) :: u, b, ls, pts
     real(dp),allocatable,dimension(:, :) :: a
     real(dp),allocatable,dimension(:, :, :) :: submat ! sub matrices for each elements, just a temporary buffer
     real(dp),allocatable,dimension(:, :)    :: subvec ! sub vectors for right hand side vector, just a temporary buffer
  end type elem
contains
  subroutine get_size(e)
    implicit none
    type(elem),intent(inout) :: e
    integer :: size

    e%nelements = 0
    e%nnodes    = 0
    size        = 0
    open(ui)
    do while(.true.)
       read(ui, *, end=10)
       size = size + 1
    end do
10  close(ui)
    e%nnodes    = size + 2
    e%nelements = size + 1
  end subroutine get_size

  subroutine allocate_arrays(e)
    implicit none
    type(elem),intent(inout) :: e
    
    allocate(e%u(e%nnodes), e%b(e%nnodes), e%pts(e%nnodes))
    allocate(e%a(e%nnodes, e%nnodes))
    allocate(e%ls(e%nelements))
    allocate(e%submat(2, 2, e%nelements))
    allocate(e%subvec(2, e%nelements))

    ! zero clear
    e%u = 0.0d0
    e%b = 0.0d0
    e%pts = 0.0d0
    e%a = 0.0d0
    e%ls = 0.0d0
    e%submat = 0.0d0
    e%subvec = 0.0d0
  end subroutine allocate_arrays

  subroutine deallocate_arrays(e)
    implicit none
    type(elem),intent(inout) :: e

    deallocate(e%u, e%b, e%a, e%ls, e%pts, e%submat, e%subvec)
  end subroutine deallocate_arrays

  subroutine read_input(e)
    implicit none
    type(elem),intent(inout) :: e
    integer :: i, size

    e%pts(1) = 0.0d0
    do i = 2, e%nnodes - 1
       read(ui, *) e%pts(i)
       if (e%pts(i) >= length) then
          write(6, *) "input value is greater than length,", length
          stop
       else if (e%pts(i) <= 0.0d0) then
          write(6, *) "input value is less than 0"
          stop
       end if
    end do
    e%pts(e%nnodes) = length
  end subroutine read_input

  subroutine calc_lengths(e)
    implicit none
    type(elem),intent(inout) :: e
    integer :: i

    do i = 1, e%nelements
       e%ls(i) = e%pts(i+1) - e%pts(i)
    end do
  end subroutine calc_lengths

  subroutine calc_submats(e)
    implicit none
    type(elem),intent(inout) :: e
    integer :: i, j, k

    ! first order element, ((x_{n+1}-x)/l_n, (x-x_n)/l_n), l_n = x_{n+1}-x_n
    ! calculated analytic derivative, not general way...
    ! need to do this by numerical derivative methods
    do k = 1, e%nelements
       do j = 1, 2
          do i = 1, 2
             if (i == j) then
                e%submat(i, j, k) =  1.0d0/e%ls(k)
             else
                e%submat(i, j, k) = -1.0d0/e%ls(k)
             end if
          end do
       end do
    end do

  end subroutine calc_submats

  ! linear function f(x) = alpha*x + beta
  function lin(alpha, beta, x) result(res)
    implicit none
    real(dp),intent(in) :: alpha, beta, x
    real(dp) :: res

    res = alpha*x + beta
  end function lin
  
  ! numerical integration
  ! gaussian quadrature for linear function
  ! integrate alpha*x + beta from a to b
  function gauss_quad(alpha, beta, a, b) result(res)
    implicit none
    real(dp),intent(in) :: alpha, beta, a, b
    real(dp) :: res, x
    integer,parameter :: npts = 2
    real(dp),dimension(npts),parameter :: weights = (/1.0d0, 1.0d0/)
    real(dp),dimension(npts),parameter :: gauss_nodes = (/-1.0d0/sqrt(3.0d0), 1.0d0/sqrt(3.0d0)/)
    integer :: i

    res = 0.0d0

    do i = 1, npts
       x   = (b - a)/2.0d0*gauss_nodes(i) + (b + a)/2.0d0
       res = res + (b - a)/2.0d0*weights(i)*lin(alpha, beta, x)
    end do
    
  end function gauss_quad

  ! calculate sub vectors of right hand side vector
  ! lambda*(\int_{x_i}^x_{i+1}dx(x_{i+1}-x)/l_i, \int_{x_i}^x_{i+1}dx(x-x_i)/l_i)
  subroutine calc_subvecs(e)
    implicit none
    type(elem),intent(inout) :: e
    integer :: i

    do i = 1, e%nelements
       e%subvec(1, i) = -1.0d0*lambda*gauss_quad(-1.0d0/e%ls(i),  1.0d0*e%pts(i+1)/e%ls(i), e%pts(i), e%pts(i+1))
       e%subvec(2, i) = -1.0d0*lambda*gauss_quad( 1.0d0/e%ls(i), -1.0d0*e%pts(i)/e%ls(i),   e%pts(i), e%pts(i+1))
    end do
    
  end subroutine calc_subvecs

  ! assemble local matrices to global one
  subroutine asm_submats(e)
    implicit none
    type(elem),intent(inout) :: e
    integer :: i, j, k

    do k = 1, e%nelements
       do j = 1, 2
          do i = 1, 2
             e%a(k+i-1, k+j-1) = e%a(k+i-1, k+j-1) + e%submat(i, j, k)
          end do
       end do
    end do

    ! the above matrix does not have an inverse matrix, so add a condition by the boundary condition: u_1 = 0
    ! | 1 0 0 0 | |u_1|   | 0 |
    ! |         | |   | = |   |
    ! |         | |   |   |   |
    ! |         | |   |   |   |
    e%a(1, 1) = 1.0d0
    e%a(1, 2) = 0.0d0
    
  end subroutine asm_submats

  ! assemble right hand side vector b
  subroutine asm_subvecs(e)
    implicit none
    type(elem),intent(inout) :: e
    integer :: i, j

    do j = 1, e%nelements
       do i = 1, 2
          e%b(j+i-1) = e%b(j+i-1) + e%subvec(i, j)
       end do
    end do

    ! u_1 = 0
    e%b(1) = 0.0d0

  end subroutine asm_subvecs

    subroutine debug(e)
    implicit none
    type(elem),intent(in) :: e
    integer :: i, j

    write(6, *) "points:"
    do i = 1, e%nnodes
       write(6, '(1pe14.5)', advance='no') e%pts(i)
    end do
    write(6, *)

    write(6, *) "lhs A:"
    do i = 1, e%nnodes
       do j = 1, e%nnodes
          write(6, '(1pe14.5)', advance='no') e%a(i, j)
       end do
       write(6, *)
    end do

    write(6, *)
    write(6, *) "rhs b:"
    do i = 1, e%nnodes
       write(6, '(1pe14.5)') e%b(i)
    end do
    
  end subroutine debug

  function exact_sol(e, x) result(res)
    implicit none
    type(elem),intent(in) :: e
    real(dp),intent(in) :: x
    real(dp) :: res

    res = lambda/2.0d0*x*(x - 2.0d0*length)
    
  end function exact_sol
  
  subroutine check(e)
    implicit none
    type(elem),intent(in) :: e
    integer :: i
    real(dp) :: exact, diff

    do i = 1, e%nnodes
       exact = exact_sol(e, e%pts(i))
       diff = e%u(i) - exact
       write(6, '(a, i, 3(1pe14.5))') "i, exact solution, 1st element FEM, diff:", i, exact, e%u(i), diff
    end do
    
  end subroutine check

  ! ax = A*x
  subroutine a_dot_x(size, a, x, ax)
    implicit none
    integer,intent(in) :: size
    real(dp),dimension(size, size),intent(in) :: a
    real(dp),dimension(size),intent(in) :: x
    real(dp),dimension(size),intent(out) :: ax
    integer :: i, j

    !$omp parallel
    !$omp workshare
    ax = 0.0d0
    !$omp end workshare
    !$omp do
    do i = 1, size
       do j = 1, size
          ax(i) = ax(i) + a(i, j)*x(j) 
       end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine a_dot_x

  ! inner product xy = x*y
  subroutine x_dot_y(size, x, y, xy)
    implicit none
    integer,intent(in) :: size
    real(dp),dimension(size),intent(in) :: x, y
    real(dp),intent(out) :: xy
    integer :: i
    
    xy = 0.0d0
    !$omp parallel do reduction(+:xy)
    do i = 1, size
       xy = xy + x(i)*y(i)
    end do
    
  end subroutine x_dot_y

  ! calc r = b - A*x
  subroutine b_minus_ax(size, a, x, b, r)
    implicit none
    integer,intent(in) :: size
    real(dp),dimension(size, size),intent(in) :: a
    real(dp),dimension(size),intent(in) :: x, b
    real(dp),dimension(size),intent(out) :: r
    real(dp),dimension(size) :: ax
    integer :: i

    call a_dot_x(size, a, x, ax)
    !$omp parallel do
    do i = 1, size
       r(i) = b(i) - ax(i)
    end do
  end subroutine b_minus_ax

  ! need to solve asymmetric matrix
  ! http://www.jicfus.jp/wiki/index.php?Bi-CGSTAB%20%E6%B3%95
  ! https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
  subroutine bicgstab(size, a, x, b)
    implicit none
    integer,intent(in) :: size
    real(dp),dimension(size, size),intent(in) :: a
    real(dp),dimension(size),intent(out) :: x
    real(dp),dimension(size),intent(in) :: b
    integer :: i, j, iter
    real(dp),dimension(size) :: r, r_new, v, v_new, p, p_new
    real(dp),dimension(size) :: r0, h, s, t
    real(dp),dimension(size) :: rr
    real(dp),dimension(size) :: xx, xx_new
    real(dp) :: rho, rho_new, omega, omega_new
    real(dp) :: alpha, beta
    real(dp) :: r0r, t2, r0v, ts, res

    !$omp parallel
    !$omp workshare
    xx = b ! initial guess
    !$omp end workshare
    !$omp end parallel
    call b_minus_ax(size, a, xx, b, r)
    !$omp parallel
    !$omp workshare
    r0 = r
    !$omp end workshare
    !$omp end parallel
    call x_dot_y(size, r0, r, r0r)

    if (r0r == 0.0d0) then
       write(6, *) "r*r0 is zero."
       stop
    end if
    
    rho   = 1.0d0
    alpha = 1.0d0
    omega = 1.0d0

    !$omp parallel
    !$omp workshare
    v = 0.0d0
    p = 0.0d0
    !$omp end workshare
    !$omp end parallel

    do iter = 1, iter_max
       res = 0.0d0
       
       call x_dot_y(size, r0, r, rho_new)
       beta = (rho_new/rho)*(alpha/omega)
       !$omp parallel do
       do i = 1, size
          p_new(i) = r(i) + beta*(p(i) - omega*v(i))
       end do
       call a_dot_x(size, a, p_new, v_new)
       call x_dot_y(size, r0, v_new, r0v)
       alpha = rho_new/r0v
       !$omp parallel do
       do i = 1, size
          h(i) = xx(i) + alpha*p_new(i)
       end do
       call b_minus_ax(size, a, h, b, rr)
       call x_dot_y(size, rr, rr, res)
       res = sqrt(res)
       if (res <= tol) then
          !$omp parallel
          !$omp workshare
          xx_new = h
          !$omp end workshare
          !$omp end parallel
          exit
       end if
       !$omp parallel do
       do i = 1, size
          s(i) = r(i) - alpha*v_new(i)
       end do
       call a_dot_x(size, a, s, t)
       call x_dot_y(size, t, s, ts)
       call x_dot_y(size, t, t, t2)
       omega_new = ts/t2
       !$omp parallel do
       do i = 1, size
          xx_new(i) = h(i) + omega_new*s(i)
       end do
       call b_minus_ax(size, a, xx_new, b, rr)
       call x_dot_y(size, rr, rr, res)
       res = sqrt(res)
       if (res <= tol) exit
       !$omp parallel do
       do i = 1, size
          r_new(i) = s(i) - omega_new*t(i)
       end do

       !$omp parallel
       !$omp workshare
       xx    = xx_new
       p     = p_new
       r     = r_new
       v     = v_new
       !$omp end workshare
       !$omp end parallel
       rho   = rho_new
       omega = omega_new

#ifdef _DEBUG
       if (mod(iter, 100) == 0) then
          write(6,*) "iter:", iter, "res:", res
          if (res.ne.res) then ! NaN check
             write(6, *) "res is NaN, iter:", iter
             stop
          end if
       end if
#endif
    end do ! iter

    if (iter>=iter_max .and. res>tol) then
       write(6, *) "did not converge."
       write(6, *) "iter_max, res:", iter_max, res
       stop
    end if

    ! converged
    write(6, *) "BiCGSTAB method converged."
    write(6, *) "iter, res:", iter, res
    !$omp parallel
    !$omp workshare
    x = xx_new
    !$omp end workshare
    !$omp end parallel
  end subroutine bicgstab
end module vars

program main
  use vars
  implicit none
  integer :: size
  integer :: i
  type(elem) :: e
  
  call get_size(e)
  write(6, *) "e%nelements:", e%nelements
  write(6, *) "e%nnodes:", e%nnodes
  call allocate_arrays(e)
  call read_input(e)
  call calc_lengths(e)
  call calc_submats(e)
  call calc_subvecs(e)
  call asm_submats(e)
  call asm_subvecs(e)
  
#ifdef _DEBUG
  call debug(e)
#endif
  ! solve Au = b
  call bicgstab(e%nnodes, e%a, e%u, e%b)
  call check(e)
  call deallocate_arrays(e)
  stop
end program main
