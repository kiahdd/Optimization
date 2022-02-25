function MultivariableNewton(~)
clc;
clear;


%% calculating roots of the gradient matrix by solving the nonlinear system of equations.
fun = @root2d;
x0 = [1,1];
x = fsolve(fun,x0)  % x1 = 2 and x2 = -1 are the roots
%% define variables 
clc;
clear;
syms x1; %% Define Variable 1
syms x2;  %% Define Variable 2

f = sqrt(x1^2+1) + sqrt(x2 ^2 +1)
%This can be extended to any number of variables, just the sizes of the matrices will be different. 
fval1=diff(f,x1) %% Define Expression 1
fval2= diff(f,x2)  %% Define Expression 2

fp=[fval1;      %% Insert them in a column matrix        
    fval2];
fp

Jac=[diff(fval1,x1),diff(fval1,x2);
     diff(fval2,x1),diff(fval2,x2)]; %% Define the jacobian in the same order
Jac    
JacT = transpose(Jac)
invert = (inv(JacT))

p1=-invert*fp  %% Term following X_n in equation 2
disp (p1)

%%
syms tol;    %% Define the tolerance vector
syms err;    %% Define the error matix that checks how closely is the result of an iteration converging on the solution.
X02 = [1;  %% Define guess values for the variables
     1];
iter=100;       %% Define Number of Iterations. More the number of iterations, longer will the code take to converge on a solution.    
tol=[1e-4;
    1e-4];     %% Specify desired tolerances.

X2=X02;        
Xold2=X02;

for i=1:iter  % 10 iterations
    h12=double(subs(p1(1),{x1,x2},{X2(1),X2(2)})); %% Evaluation of the partial differentials with the guess values.
    h22=double(subs(p1(2),{x1,x2},{X2(1),X2(2)}));
    d=[h12;
       h22];
    d = double(d)
    X2=X2+d
    subs( f, {x1,x2},{X2(1),X2(2)})
    err=double(abs(X2-Xold2))
    Xold2=X2;
    if (err(1)<tol(1) && err(2)<tol(2))
        break
    end
    disp("******")
    disp(i)
end

end
