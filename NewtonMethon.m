function NewtonMethon(~)
% x1 = -10;
% x2 = 10;

% fun =  @(x,y)(x-2).^4+(x-2).^2*y^2+(y+1).^2;
% min = fminbnd(fun,x1,x2)
% x = linspace(-5,5)
% plot(x, fun(x))
clear; clc;
f =  @(x)x.^4-32.*x.^2;
fp = @(x) 4.*x.^3-64.*x;
fpp = @(x) 12.*x.^2-64;

x0 = 2;

N = 20; 
tol = 1E-10;

x(1) = x0; % Set initial guess
n = 2; 
nfinal = N + 1; 

while (n <= nfinal)
  fe = fp(x(n - 1));
  fpe = fpp(x(n - 1));
  x(n) = x(n - 1) - fe/fpe;
  if (abs(fe) <= tol)
    nfinal = n; 
    break;
  end
  n = n + 1;
end
    disp(x)
    fp(x)
    f(x)
end