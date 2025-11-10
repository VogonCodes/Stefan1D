# 1D Stefan problem
The analytical solution in the solidified layer for dimensionless temperature

$$ \Theta = \frac{T-T_{\mathrm{top}}}{T_{\mathrm{melt}}-T_{\mathrm{top}}} = \frac{\mathrm{erf}\eta}{\mathrm{erf}\lambda},$$

where $\eta=\frac{y}{2\sqrt{\kappa t}}$ is the dimensionless coordinate and $\lambda=\eta_{m}=\frac{y_{m}}{2\sqrt{\kappa t}}$ is the dimensionless position of the solidification boundary, which is constant over time.

The speed of solidification boundary is

$$ \frac{\mathrm{d}y_{m}}{\mathrm{d} t} = \lambda\sqrt{\frac{\kappa}{t}}, $$

but from the conservation of energy also

$$ \rho L\frac{\mathrm{d} y_{m}}{\mathrm{d} t} = k\frac{\partial T}{\partial y}\Bigg|_{y=y_{m}}. $$

Using the analytical solution for temperature, one can get to the transcendental equation for $\lambda$ (should be constant over time!)

$$ \frac{L\sqrt{\pi}}{c(T_{\mathrm{melt}}-T_{\mathrm{top}})} = \frac{\exp(-\lambda^{2})}{\lambda\mathrm{erf}\lambda} $$

## Numerical solution
Solve the heat equation

$$ \frac{\partial T}{\partial t} = \kappa\frac{\partial^2 T}{\partial x^2} $$

with boundary conditions

$$
T(0,t) = T_{\mathrm{top}}
\,,\qquad
T(y_{m},t) = T_{\mathrm{melt}}
\,,
$$

where the change of $y_{m}$ is described by

$$ \rho L\frac{\mathrm{d} y_{m}}{\mathrm{d} t} = k\frac{\partial T}{\partial y}\Bigg|_{y=y_{m}}. $$

### Heat equation
#### Grid
The grid is equidistant.
After calculating $\Delta y_{m}$, recaltulate grid nodes position and $\Delta x$.
If the new $\Delta x>\Delta x_{\mathrm max}$, change number of nodes by 1 (i.e. $N$ -> $N+1$).
Not the best way, but it works fine for me (but it wouldn't be so hard to break it).
Then calculate temperature in new nodes from the old grid by linear interpolation.

#### Time step
Use CFL criterion

$$ \Delta t=\mathrm{CFL}\frac{(\Delta x)^{2}}{\kappa} $$

#### Central differences
At time $t_{n+1}$ the solution is

$$ T^{n+1}\_{j} = T^{n}\_{i} + \mu(T^{n}\_{j+1} - 2T^{n}\_{j} + T^{n}_{j-1})
\qquad j\in\{1,\dots,N-1\}\,, $$

where $\mu=\frac{\kappa\Delta t}{(\Delta x)^{2}}$,
with $T^{n}\_{0}=T_{\mathrm{top}}$, $T^{n}\_{J}=T_{\mathrm{melt}}$ $\forall n\in\{0,\dots,N\}$.

I just put that into a `for` loop and I don't regret it.

#### Crank-Nicholson

$$ -\mu T^{n+1}_{j+1} + 2(1+\mu) T^{n+1}_{j} - \mu T^{n+1}_{j-1} = \mu T^{n}_{j+1} + 2(1-\mu)T^{n}_{j} + \mu T^{n}_{j-1} \qquad j\in\{1,\dots,J-1\} $$

with $T^{n}\_{0}=T_{\mathrm{top}}$, $T^{n}\_{J}=T_{\mathrm{melt}}$ $\forall n\in\{0,\dots,N\}$.

This is a system of $J-2$ linear equations

$$
\begin{pmatrix}
    2(1+\mu)    & -\mu      & 0         & 0       & 0 & \dots & 0       & 0         \\
    -\mu        & 2(1+\mu)  & -\mu      & 0       & 0 & \dots & 0       & 0         \\
    0           &-\mu       & 2(1+\mu)  & -\mu    & 0 & \dots & 0       & 0         \\
    \vdots      &           &           & \ddots  &   &       & 0       & 0         \\
    0           & \dots     & \dots     &         &   &       & -\mu    & 2(1+\mu)  
\end{pmatrix}
\begin{pmatrix}
    T^{n+1}_{1}\\
    T^{n+1}_{2}\\
    T^{n+1}_{3}\\
    \vdots\\
    T^{n+1}_{J-1}
\end{pmatrix}=
\begin{pmatrix}
\mu T^{n}_{2} + 2(1-\mu)T^{n}_{1} + 2\mu T^{n}_{0}\\
\mu T^{n}_{3} + 2(1-\mu)T^{n}_{2} + \mu T^{n}_{1}\\
\mu T^{n}_{4} + 2(1-\mu)T^{n}_{3} + \mu T^{n}_{2}\\
\vdots\\
2\mu T^{n}_{J} + 2(1-\mu)T^{n}_{J-1} + \mu T^{n}_{J-2}
\end{pmatrix}
$$

##### Thomas algorithm
Basically Gauss for tridiagonal matrices.
If $a_{i}$ are subdiagonal elements, $b_{i}$ diagonal and $c_{i}$ are superdiagonal, then individual equations have the form of

$$
a_{j}T^{n+1}_{j+1}+b_{j}T^{n+1}_{j}+c_{j}T^{n+1}_{j-1}=d_{j}.
$$

First, the matrix is transformed to upper triangular with unit diagonal by Gauss: 

$$
a_{i}'=0
\qquad
b_{i}'=1
\qquad
c_{i}=\frac{c_{i}}{b_{i}-a_{i}c_{i-1}'}
\qquad
d_{i}=\frac{d_{i}-a_{i}d_{i-1}'}{b_{i}-a_{i}c_{i-1}'},
$$

meaning the equations are in the form of

$$ T^{n+1}_{j}+c_{j}' T^{n+1}_{j-1}=d_{j}'\,, $$

and then using backsubstitution (start from either $j=1$ or $j=J-1$, where $T_{j}^{n+1}=d_{j}'$ and then use the equation above.

### Solidification boundary
Use forward Euler for $y_{m}$ and backward Euler for $\nabla T$

$$ (y_{m})^{n+1}=(y_{m})^{n} + \frac{k\Delta t}{\rho L}\frac{T^{n}_{J}-T^{n}_{J-1}}{\Delta x} $$

## Analytical solution
### Solidification boundary
Newton to find $\lambda$:

$$
\begin{aligned}
    f &= \frac{\exp(-\lambda^{2})}{\lambda\mathrm{erf}\lambda} - \frac{L\sqrt\pi}{c\Delta T}\\
    f' &= -\frac{\exp(-\lambda^{2})}{\lambda^{2}\mathrm{erf}\lambda}-\frac{2\exp(-\lambda^{2})}{\sqrt{\pi}\lambda\mathrm{erf}^{2}\lambda}
\end{aligned}
$$

Thus,

$$ \lambda=\lambda-\frac{f}{f'} $$

### Temperature
If $y<y_{m}$, the temperature is given by

$$ T_{j}^{n} = T_{\mathrm{top}} + \Delta T \frac{\mathrm{erf}{\eta_{j}^{n}}}{\mathrm{erf}\lambda}, $$

where $\eta_{j}^{n}=\frac{j\Delta x}{2\sqrt{\kappa t_{n}}}$


---

jinak:
- resit na gridu vetsi, nez led -> v nodu tesne pod solidifikacni frontou natvrdo predepsat Tmelt (funguje at least pro konecne diference, pro CN by imho melo taky)
