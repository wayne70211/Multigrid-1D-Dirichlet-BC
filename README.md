# 1D Poisson Multigrid Solver

### 1D Poisson Equation Problem
<br>
<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\large&space;\nabla^2&space;u&space;=&space;f(x)&space;=&space;\frac{\sin&space;\pi&space;x&space;&plus;&space;\sin&space;16\pi&space;x}{2}" title="\large \nabla^2 u = f(x) = \frac{\sin \pi x + \sin 16\pi x}{2}" />
</p>
<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\large&space;x\in&space;[0,1]" title="\large x\in [0,1]" /> 
</p>

with *Dirichlet Boundary Condition*  <br>
<br>
<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\large&space;u(0)=0\;&space;\;&space;u(1)=1" title="\large u(0)=0\; \; u(1)=1" /><br>
</p>

The exact solution is <br>
<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\large&space;u(x)=-\frac{\sin&space;\pi&space;x}{2\pi^2}&space;-&space;\frac{\sin&space;16\pi&space;x}{512\pi^2}" title="\large u(x)=-\frac{\sin \pi x}{2\pi^2} - \frac{\sin 16\pi x}{512\pi^2}" /> <br>
</p>

---
### Exact Solver on Coarsest Level
The coarsest grid should be solved to exact solution<br>
The program provides two solvers
1. **Gauss-Seidel** with residual L2-norm less than tolerance
2. ~~**LU decomposition**~~ (Not include)
3. ~~**CG Solver**~~ (Not include)

``` fortran
tol = 1.0d-8
exact_solver = 1 
```

###  Set Level

Set arbitrary level. The coarsest grid should greater than 3

``` fortran
if (Level_Nx(Level_num).le.3) then
    write(*,*) Level_num," level is high for ",Nx," grid "
    stop
end if
```

###  Set Multigrid Solver

**Gauss-Seidel** scheme is selected to be `smoother` <br>

**V Cycle Scheme**
``` fortran
subroutine MG_Vcycle(Nx,dx,RHS,U,Level_num,tol,exact_solver,I_cycle)
    implicit none
    integer,intent(in) :: Nx,Level_num,exact_solver
    real*8 ,intent(in) :: dx,tol
    real*8,dimension(0:Nx), intent(in)     :: RHS
    real*8,dimension(0:Nx), intent(inout)  :: U
    integer,dimension(Level_num) :: Level_Nx
    real*8 ,dimension(Level_num) :: Level_dx
    real*8  :: rms0,rms,rmsc
    integer :: Max_Iter,Relax_Iter
    integer :: i,k,Level,iter
    integer, intent(inout) :: I_cycle
```

---
## Multigrid V cycle Result
The following tables are the result of **Multigrid V cycle** with 7, 8, and 9 levels <br>

|   Solver  |  Np   |  Cycle  |   CPU Time  |  
| :---:     | :---: | :---:   |   :---:     |  
| MG7       | 4096  |   7     |     0.0085  | 
| MG8       | 4096  |   7     |     0.0021  | 
| MG9       | 4096  |   7     |     0.0011  | 

|   Solver  |  Np   |  Cycle  |   CPU Time  |  
| :---:     | :---: | :---:   |   :---:     |
| MG7       | 8192  |   7     |     0.0673  | 
| MG8       | 8192  |   7     |     0.0088  | 
| MG9       | 8192  |   7     |     0.0030  | 

|   Solver  |  Np   |  Cycle  |   CPU Time  |  
| :---:     | :---: | :---:   |   :---:     |
| MG7       | 16384 |   7     |     0.5301  | 
| MG8       | 16384 |   7     |     0.0684  | 
| MG9       | 16384 |   7     |     0.0110  | 

## Compile and Run 
Use [`PGI Compiler`](https://www.pgroup.com/products/community.htm)

```shell
make
```

Run
```shell
./Run
```
## Reference
See [wikipedia](https://en.wikipedia.org/wiki/Multigrid_method)
