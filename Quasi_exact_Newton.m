function Quasi_exact_Newton(~)


clc; clear;

Q = [5 -3; -3 2]
b = [0; 1]
x0 = [0;0]
H0 = eye(2)

f = @(x) 1/2 * transpose(x) * Q * x - transpose(b) * x + log2(pi);

[x,fval,exitflag,output,grad,hessian] = fminunc(f,x0)

fp = @(x) Q * x - b;

d0 = - inv(H0) * fp(x0)
x1 = @(alpha) x0 + alpha * d0;

fp(x0)'*d0
(d0'*Q)
(d0'*Q)*d0
alpha0 = -(fp(x0)'*d0)/((d0'*Q)*d0)

x1 = x1(alpha0)
f(x1)

pause
S = x1-x0
y = fp(x1) - fp(x0)

% H1_inv = (1 - (S*y')/(y'*S))*inv(H0)*(1-(y*S')/(y'*S))+(S*S')/(y'*S)

H1 = H0 + (y*y')/(y'*S)-(H0*S*S'*H0')/(S'*H0*S)

H1_inv = inv(H1)
 disp('whether the value below is equal to S')
 
H1_inv * y
S

d1 = - H1_inv * fp(x1)

x2 = @(alpha) x1 + alpha * d1;

alpha1 = -((Q*x0-b)'*d1)/(d1'*Q*d1)

x2 = x2(alpha1)

f(x2)

end