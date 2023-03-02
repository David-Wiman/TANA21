function [x, k] = SuccessiveOverRelaxation(A,b,w,x,tol)
%  Gauss-Seidel iterative method to approximate the solution of a
%  linear system Ax=b up to a user defined tolerance
%
%  INPUT: 
%    A   - n by n square, non-singular, diagonally dominant matrix
%    b   - n by 1 right hand side vector
%    w   - 1 by 1 constant
%    x   - n by 1 vector containing that initial guess for the iteration
%    tol - user set tolerance for the stopping condition in the iteration 
%
%  OUTPUT:
%    x - n by 1 vector containing the iterative solution
%    k - number of iterations

%  get the system size
   n = length(A);
   
%  defining the maximum number of iterations   
   max_its = 22000;
   
%  looping through the iteration index k
   for k = 1:max_its
%  looping through the rows of the x vector
       for i = 1:n
           sum_1 = 0;
           sum_2 = 0;
%  creating the two sums needed in the formula           
           for j = 1:i-1
                sum_1 = sum_1 + A(i,j)*x(j,k+1);
           end
           for j = i+1:n
               sum_2 = sum_2 + A(i,j)*x(j,k);
           end
%   updating the iterative solution           
           x(i,k+1) = (1-w)*x(i,k)+(b(i)-sum_1-sum_2)*(w/A(i,i));
       end
%   checking if the tolerance i met and stopping if that's the case       
       r = b - A*x(:,k);
       if norm(r,2) <= tol*norm(b,2)
           break
       end
   end

%   isolating the last column of x which contains the iterative solution
%   that met the stopping condition
    x = x(:,end);

end