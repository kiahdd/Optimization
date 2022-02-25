function NewtonMethodBacktracking(~)
clc;
flag        = input('inital values at zero? [y]/[n]  ','s')
if flag == 'n'
    x0      = input('input new initial values:  ');
    check   = size(x0);
    if check(2) == 2
        x0  = transpose(x0);
    end
else
    x0  = [0;0];
end
%% define variables 
syms x1; %% Define Variable 1
syms x2;  %% Define Variable 2
f       = input('input function:  ');
fp      = [diff(f,x1);diff(f,x2)]

Jac     = [diff(fp(1),x1),diff(fp(1),x2);
            diff(fp(2),x2),diff(fp(2),x2)] %% Define the jacobian in the same order
 
JacT    = transpose(Jac)
invert  = (inv(JacT))

betha   = input('set the Betha Value (0,1):   ');
gamma   = input('set the Gamma Value (0,1):   ');
dd      = invert*fp  %% Term following X_n in equation 2

%%
syms tol;    %% Define the tolerance vector
syms err;    %% Define the error matix that checks how closely is the result of an iteration converging on the solution.

iter        = 100;       %% Define Number of Iterations. More the number of iterations, longer will the code take to converge on a solution.    
tol         = [1e-10;
              1e-10];     %% Specify desired tolerances.

xK          =  x0;        
xold        =  xK;
alphaMax    = 10; % this is the maximum step length
alpha_K       = alphaMax;

for K=1:iter  % 10 iterations
    disp('starting iteration: ')
    disp(K)
    dx1     = double(subs(dd(1),{x1,x2},{xK(1),xK(2)})); %% Evaluation of the partial differentials with the guess values.
    dx2     = double(subs(dd(2),{x1,x2},{xK(1),xK(2)}));
    dx      =   [dx1;dx2]; %dx is a 2 by 1 vector
    
    while true

        fold    = double(subs(f,{x1,x2},{xold(1),xold(2)}));                         % THIS IS f(x)
        fnew    = double(subs(f,{x1,x2},{xK(1)-alpha_K*dx(1),xK(2)-alpha_K*dx(2)})); % THIS IS f(x+alpha*dx)
        fp_old1 = double(subs(fp(1),{x1,x2},{xold(1),xold(2)}));                      % THIS is Jacobian f (x)
        fp_old2 = double(subs(fp(2),{x1,x2},{xold(1),xold(2)}));
        fp_old  = [fp_old1;fp_old2];
        
        
        if fold - fnew < - gamma *  alpha_K * transpose(fp_old) * dx
           alpha_K = betha * alpha_K;
           disp('reducing alpha to: ')
           disp(alpha_K)
        else 
            break;
        end
    end
    xK      = xK - alpha_K * dx
    err     = double(abs(xK-xold));
    xold    = xK;
    
    if (err(1)<tol(1) && err(2)<tol(2))
        break
    end
    disp("******")
    




end
disp('Optimal points of the function are: ')
disp(xold)
disp('final optimal step size is: ')
disp(alpha_K)