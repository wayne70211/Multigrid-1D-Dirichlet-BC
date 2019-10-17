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
