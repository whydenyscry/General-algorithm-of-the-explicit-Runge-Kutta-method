# General algorithm of the explicit Runge—Kutta method
Below I will describe an algorithm for solving IVP by any explicit Runge-Kutta method of any order for any dimensionality of the system, I have not seen such an implementation anywhere. 
The advantage of this algorithm is that the code will look much more compact, since you will not have to create cumbersome and long expressions to compute this or that expression for each individual IVP. 
All you have to do is fill in the Butcher table for the method you want the IVP to be solved by.

## Explicit Runge—Kutta methods. Butcher tableau
Let an initial value problem be specified as follows:

$$ \dot{\mathbf{x}}=\mathbf{f}\left(t,\mathbf{x}\right),\quad t \in \left[t_0,t_\text{end}\right],\quad \mathbf{x}\left(t_0\right) = \mathbf{x}_0 \in \mathbb{R}^m, $$

where $\mathbf{x}=\left[x_1,\dots,x_m\right]^\mathbf{T},\quad
	\mathbf{f}\left(t,\mathbf{x}\right)=\left[f_1\left(t,x_1,\dots,x_n\right),\dots,f_m\left(t,x_1,\dots,x_n\right)\right]^\mathbf{T}.$
	
The $s$-stage Runge-Kutta method can be expressed as follows:

$$ \mathbf{x}_{n+1} = \mathbf{x}_n + \tau \sum_{i=1}^{s} b_i \mathbf{k}_{i}^{(n)} $$
where
$$ \begin{cases}
\mathbf{k}_{1}^{(n)} = \mathbf{f}(t_n, \mathbf{x}_n), \\
\vdots \\
\mathbf{k}_{i}^{(n)} = \mathbf{f}\left(t_n + c_i \tau, \mathbf{x}_n + \tau \sum_{j=1}^{i-1} a_{i,j} \mathbf{k}_{j}^{(n)}\right), \quad i = \overline{2,s}.
\end{cases} $$

	
## References
1. Butcher, J. C. (2016). Numerical methods for ordinary differential equations. John Wiley & Sons.
