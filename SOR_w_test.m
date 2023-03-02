function w = SOR_w_test(A,b,x,tol)
%  A test function to experimentally determine the optimal value of omega
%  for the Successive Over-relaxation method
%
%  INPUT: 
%    A   - n by n square, non-singular, diagonally dominant matrix
%    b   - n by 1 right hand side vector
%    x   - n by 1 vector containing that initial guess for the iteration
%    tol - user set tolerance for the stopping condition in the iteration
%
%  OUTPUT
%    w   - the omega that results in the least number of iterations

%   allocates space for the vector containing the number of iterations for
%   each omega
    iteration_vector = zeros(1,101);
    
%   looping through the different omegas with 0.01 increments and computing
%   and storing the number of iterations needed to met the tolerance   
    for w = 1:0.01:2
        [~,k] = SuccessiveOverRelaxation(A,b,w,x,tol);
        index = round((w-1)*100+1);
        iteration_vector(index) = k;
    end
    
%   finding the minimum number of iterations needed and finding the
%   corresponding omega   
    [~,min_index] = min(iteration_vector);
    w = [1:0.01:2];
    w = w(min_index);
    
end