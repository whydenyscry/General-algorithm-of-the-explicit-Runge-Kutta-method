# General algorithm of the explicit Runge—Kutta method
Below I will describe an algorithm for solving IVP by any explicit Runge-Kutta method of any order for any dimensionality of the system, I have not seen such an implementation anywhere. 
The advantage of this algorithm is that the code will look much more compact, since you will not have to create cumbersome and long expressions to compute this or that expression for each individual IVP. 
All you have to do is fill in the Butcher table for the method you want the IVP to be solved by.

## Explicit Runge—Kutta methods. Butcher tableau

$$	\dot{\mathbf{x}}=\mathbf{f}\left(t,\mathbf{x}\right),\quad t \in \left[t_0,t_\text{end}\right],\quad \mathbf{x}\left(t_0\right) = \mathbf{x}_0 \in \mathbb{R}^m.$$
## References
1. Butcher, J. C. (2016). Numerical methods for ordinary differential equations. John Wiley & Sons.
