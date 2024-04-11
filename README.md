# General algorithm of the explicit Runge-Kutta method
Below I will describe an algorithm for solving IVP by any explicit Runge-Kutta method of any order for any dimensionality of the system, I have not seen such an implementation anywhere. 
The advantage of this algorithm is that the code will look much more compact, since you will not have to create cumbersome and long expressions to compute this or that expression for each individual IVP. 
All you have to do is fill in the Butcher table for the method you want the IVP to be solved by.
