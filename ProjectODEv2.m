function [y_approx, error] = ProjectODEv2(num_nodes, tol)
% A Forward Euler method for solving the differential equation
% y'' = cos(x)y'-sin(x)y on the interval [0,3*pi/2] with the initial
% conditions y(0)=1 and y(3*pi/2)=exp(-1)
% Consider running the SOR_w_test.m after the first run to get the optimal
% w for your matrices
%
%   INPUT
%      num_nodes  - the number of nodes used for the approximation
%      tol        - user set tolerance for the stopping condition in the
%                   iteration
%
%   OUTPUT
%      y_approx   - n by 1 vector with the approximate function values to y
%      error      - the error between the exact and approximate function
%                   values calculated in the infinity norm

%   adds the Matrix Solver directory to the current frame (only on my computer)
    addpath 'C:\Users\David Wiman\OneDrive\Dokument\MATLAB\TANA21\Matrix Solver'

%   sets the interval, time step length, size of the system and allocates
%   space for A
    h = (3*pi)/(2*num_nodes);
    x = [0:h:(3*pi)/2];
    n = length(x);
    A = zeros(n,n);
    
%   computes the exact function values
    y_exact = exp(sin(x))';
    
%   computes and inserts the correct elements into A on row 2 to n-1
    for i = 2:n-1
        alfa = (1/h^2)+(cos(i*h)/(2*h));
        beta = (-2/(h^2))+sin(i*h);
        gamma = (1/h^2)-(cos(i*h)/(2*h));
        A(i,i-1:i+1) = [alfa beta gamma];
    end

%   computes and inserts the correct elements into A for row 1 and n
    A(1, 1:3) = [(-5)/(h^2)-(2*cos(h))/h (4/(h^2))+(cos(h))/(2*h) (-1/(h^2))];
    A(n,n-2:n) = [(-1)/(h^2) (4/(h^2))-(cos(n*h)/(2*h)) (-5)/(h^2)+(2*cos(n*h))/h];

%   creats the right hand side vector and inserts the correct elements
    b = zeros(n,1);
    b(1,1) = (-2)/(h^2)-(3*cos(h))/(2*h)-sin(h);
    b(n,1) = ((-2)/(h^2)+(3*cos(n*h)/(2*h))-sin(n*h))*exp(-1);

%   solves the maxtrix system iteratively using the SOR method
    [y_approx,~] = SuccessiveOverRelaxation(A,b,1.95,ones(n,1),tol);

%   computes the error between the approximate and the exact solution in
%   the infinity norm
    error = norm(y_exact-y_approx, inf);
    
%   plots the exact and approximate solutions in the same graph to visually
%   confirm their resemblens
    plot(x,y_exact,'b');
    hold;
    plot(x,y_approx,'r--');
    title("Exact & Approximate ODE Solutions for X Nodes") 
    xlabel("t")
    ylabel("y")
    legend("Exact solution","Approximate solution")
    
end