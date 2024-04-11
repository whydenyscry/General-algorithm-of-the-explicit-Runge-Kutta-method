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

$$ \mathbf{x}_{n+1} = \mathbf{x}_n+\tau\sum\limits\_{i=1}^{s}b\_i\mathbf{k}\_{i}^{(n)}, $$

where 

$$ 
\begin{cases}
			\mathbf{k}\_{1}^{(n)} = \mathbf{f}\left(t_n,\mathbf{x}\_n\right),\\
			\vdots\\
			\mathbf{k}\_{i}^{(n)} = \mathbf{f}\left(t_n + c_i \tau, \mathbf{x}\_n + \tau\displaystyle\sum_{j=1}^{i-1} a_{i,j}\mathbf{k}_{j}^{(n)}\right),
\end{cases}
$$

$$ i=\overline{2,s}, $$

$\tau$ — time discretization step.

Method coefficients are conveniently set in the form of a **Butcher tableau**:

$$ 
\begin{array}{r|c}
			\mathbf{c} & \mathbf{A} \\
			\hline
			& \mathbf{b}^{\mathbf{T}}
		\end{array} \quad \Rightarrow
		\begin{array}{r|ccccc}
			0     &         &         &         & \\
			c_2   & a_{2,1}  &         &         & \\
			c_3   & a_{3,1}  & a_{3,2}  &         & \\
			\vdots& \vdots  & \vdots  & \ddots  & \\
			c_s   & a_{s,1}  & a_{s,2}  & \cdots  & a_{s,s-1} \\
			\hline
			& b_1     & b_2     & \cdots  & b_{s-1} & b_s
		\end{array},
$$

$$ \mathbf{c},\mathbf{b} \in \mathbb{R}^s,\quad \mathbf{A} \in \mathbb{R}^{s\times s}. $$

In the program implementation other elements of the matrix $\mathbf{A}$ are given by zeros, for example, the Butcher table for the classical method of order 4 is given in the program as follows:

$$ \begin{array}{r|cccc}
			0          &             &             &           & \\
			1/2     & 1/2        &             &           & \\
			1/2     & 0          & 1/2        &           & \\
			1           & 0          & 0           & 1        & \\
			\hline
			& 1/6      & 1/3      & 1/3      & 1/6    \\
\end{array} \quad \Rightarrow \quad \begin{array}{c|cccc}
				0   & 0   & 0   & 0   & 0\\
				1/2 & 1/2 & 0   & 0   & 0\\
				1/2 & 0   & 1/2 & 0   & 0\\
				1   & 0   & 0   & 1   & 0\\
				\hline
				& 1/6 & 1/3 & 1/3 & 1/6\\
		\end{array} $$
		
## Description of the implemented algorithm
Of course, you can implement the algorithm described in the last section as well, and it will work the same way as the algorithm I will describe below. 

So, the algorithm is based on the application of general matrix algebra:

$$ \mathbf{x}_{n+1} = \mathbf{x}_n+\tau\mathbf{K}^{(n)}\mathbf{b}.$$

To begin with, at each iteration we need to initialize the matrix $\mathbf{K}^{(n)}$ of the corresponding size as a zero matrix and this matrix is interpreted as follows:

$$ \mathbf{K}^{(n)}\_{m\times s}=\left[\mathbf{k}_1^{(n)},\mathbf{k}_2^{(n)},\ldots,\mathbf{k}_s^{(n)}\right]=\mathbf{0}\_{m\times s}, $$

and the matrix $\mathbf{A}$:

$$ \mathbf{A}\_{s\times s} = 
		\begin{bmatrix}
			\mathbf{a}_1^\mathbf{T}
			\\
			\mathbf{a}_2^\mathbf{T}
			\\
			\vdots 
			\\
			\mathbf{a}_s^\mathbf{T}
		\end{bmatrix}.
$$

Then the formulas for filling the matrix $\mathbf{K}^{(n)}$ can be represented as follows:
$$ \begin{cases}
			\mathbf{k}\_{1}^{(n)} = \mathbf{f}\left(t_n,\mathbf{x}_n\right),\\
			\vdots\\
			\mathbf{k}\_{i}^{(n)} = \mathbf{f}\left(t_n + c_i \tau, \mathbf{x}_n + \tau\mathbf{K}^{(n)}\_{m\times i-1}\mathbf{a}_i\right),
	\end{cases}
$$

$$ i=\overline{2,s}.$$



## References
1. Butcher, J. C. (2016). Numerical methods for ordinary differential equations. John Wiley & Sons.
