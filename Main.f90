!-----------------------------------------------------------------------------!
!  Multigrid method for solving 1D Poisson equation
!     d2u/dx2 + d2u/dy2 = f(x,y)
!     Drichlet b.c.
!-----------------------------------------------------------------------------!

program Poisson_1D
  use omp_lib
  use Multigrid_1D_Dirichlet_BC
  use Multigrid_1D_Dirichlet_BC_OMP
  implicit none
  integer::i,k,nx,exact_solver,Level_num,I_cycle
  character(len=7) :: mode
  character(len=1) :: Level,version
  real*8,dimension(:),allocatable ::u,f,ue,e
  real*8,dimension(:),allocatable ::x
  real*8 ::dx,pi,tol,rms,x0,xL,start,finish

  !Domain
  x0 = 0.0d0 !left
  xL = 1.0d0 !right

  pi = 4.0d0*datan(1.0d0)

  !Tolerance
  tol = 1.0d-8

  !---------------------------------------------------------------------------!
  ! Exact Solver
  !     1. Gauss Seidel Solver
  !     2. LU decomposition [Not support]
  !     3. CG Solver        [Not support]
  !---------------------------------------------------------------------------!
  exact_solver = 1  !
  !---------------------------------------------------------------------------!

  open(66,file='Result.txt')
  write(66,*)"-------------------------------------------------------"
  write(66,*)"Solver     Np     Iter        CPU Time        L2-norm"
  write(66,*)"-------------------------------------------------------"

  write(*,*)"-------------------------------------------------------"
  write(*,*)"Solver     Np     Iter        CPU Time        L2-norm"
  write(*,*)"-------------------------------------------------------"


  do k=12,14

    ! Multigrid Level
    do Level_num=7,9

    !number of points
    nx = 2**k !number of grid points in x

    !grid spacing (spatial)
    dx = (xL-x0)/dfloat(nx)

    !spatial coordinates
    allocate(x(0:nx))
    do i=0,nx
    x(i) = x0 + dfloat(i)*dx
    end do

    allocate(u(0:nx))
    allocate(f(0:nx))
    allocate(e(0:nx))
    allocate(ue(0:nx))

    !---------------------------------------------!
    ! Exact solution
    !---------------------------------------------!
    do i=0,nx
    f(i) = 0.5d0*(dsin(pi*x(i))+dsin(16.0d0*pi*x(i)))
    ue(i)= -dsin(pi*x(i))/(2.0d0*pi*pi)-dsin(16.0d0*pi*x(i))/(512.0d0*pi*pi)
    end do


    !Numerical solution:
    do i=0,nx
    u(i)=0.0d0
    end do

    !Boundary conditions
    u(0)  = 0.0d0
    u(nx) = 1.0d0

    !---------------------------------------------------------------------------!
    ! Solver:
    !---------------------------------------------------------------------------!

    !Level_num = 6
    write(Level,'(I1.1)') Level_num
    write(version,'(I1.1)') exact_solver
    mode = 'MG'//Level//'-V'//version

    ! Numerical solution:
    do i=0,nx
    u(i)=0.0d0
    end do

    start=omp_get_wtime()
    call MG_Vcycle(Nx,dx,F,U,Level_num,tol,exact_solver,I_cycle)
    finish=omp_get_wtime()

    !----------------------!
    ! Error analysis:
    !----------------------!
    do i=0,nx
    e(i) = dabs(u(i)-ue(i))
    end do

    !L-2 Norm:
    call L2norm(nx,e,rms)

    write(*,'(A7,i8,i8,f18.12,f15.10)') mode,nx,I_cycle,finish-start,rms
    write(66,'(A7,i8,i8,f18.12,f15.10)') mode,nx,I_cycle,finish-start,rms

    !---------------------------------------------------------------------------!
    mode = 'PMG'//Level//'-V'//version

    ! Numerical solution:
    do i=0,nx
    u(i)=0.0d0
    end do

    start=omp_get_wtime()
    call P_MG_Vcycle(Nx,dx,F,U,Level_num,tol,exact_solver,I_cycle)
    finish=omp_get_wtime()

    !----------------------!
    ! Error analysis:
    !----------------------!
    do i=0,nx
    e(i) = dabs(u(i)-ue(i))
    end do

    call L2norm(nx,e,rms)

    write(*,'(A7,i8,i8,f18.12,f15.10)') mode,nx,I_cycle,finish-start,rms
    write(66,'(A7,i8,i8,f18.12,f15.10)') mode,nx,I_cycle,finish-start,rms

    !---------------------------------------------------------------------------!
    write(*,*)"-------------------------------------------------------"
    write(66,*)"-------------------------------------------------------"

    deallocate(u,f,e,ue,x)
    enddo

  enddo

  close(66)

end program
