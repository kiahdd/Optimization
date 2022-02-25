function gradientDecsentBackTracking()
% this function applies gradient descent method with backtracking
% to a costum function f that needs to be defined by the user
% the user defined function must have only two variables e.g x1, x2
% example enter function: >> x1^2 + x2^2
% initial values for x are defined by the user
% example enter intial points: >> [2;1]



clc; clear; close all;
% Define Variable 1, Variable 2
syms x1 x2;

% ask the user to input the function and set initial conditions

fprintf('systems variables are %s and %s \n', x1,x2)
f           = input('input function e.g x1^2+x2^2:      ');
X0         = input('input initial points e.g [2;1]:      ');
% creating the name of the outful excel file
fileName = input('input output file name  e.g: test.xlsx:     ','s');
path = fullfile('.','tables');

if ~exist(path, 'dir')
    mkdir(path)
end
FileName = fullfile(path,fileName)

% define Jacobian matrix for applying gradient descent method
Jac        = [diff(f,x1);
    diff(f,x2)];
disp(Jac)

% default values of gradient descent method parameters
betha      = 0.5;
gamma    = 0.25;
s            = 2;         % this is the maximum step length
tol          = 1e-5;    % tolerence
iter         = 1000;    % maximum number of iterations

Xk          =  X0;      % starting value of x at iteration k
Xold        =  Xk;      % value of x at iteration k-1

% descent direction of gradient descent method
descentDirection = -Jac;

% initialize alpha with s
alpha_K              = s;

% initialize function value
functionValue      = 1e+5;


% Table definition to store variables
t       = array2table(zeros(0,7));

t.Properties.VariableNames = {'iteration','Initial Point', ...
    'Search Direction', 'Step length','new point','error',...
    'function value'};


for K=1:iter  % 10 iterations
    fprintf ('starting iteration: %d \n', K)
    
    
    % Evaluation the descent direction at x_K
    DescentDirection_Xk  = ...
        [ double(subs(descentDirection(1),{x1,x2},{Xk(1),Xk(2)}));
        double(subs(descentDirection(2),{x1,x2},{Xk(1),Xk(2)}))];
    
    % Backtracking algorithm
    while true
        
        % THIS IS f(x)
        fold        = double(subs(f,{x1,x2},{Xold(1),Xold(2)}));
        
        % THIS IS f(x+alpha*dx)
        fnew       = double(subs(f,{x1,x2}, ...
            {Xk(1) + alpha_K * DescentDirection_Xk(1),...
            Xk(2) + alpha_K * DescentDirection_Xk(2) ...
            }));
        
        % THIS is Jacobian f (x)
        Jac_old =  [ double(subs(Jac(1),{x1,x2},{Xold(1),Xold(2)}));
            double(subs(Jac(2),{x1,x2},{Xold(1),Xold(2)})) ];
        
        
        % check the backtracking condition
        if fold - fnew < - gamma *  alpha_K * transpose(Jac_old) * DescentDirection_Xk
            % reduce alpha
            alpha_K = betha * alpha_K;
            fprintf('reducing alpha to %f \n', alpha_K)
            
        else
            % break the backtracking algorithm
            functionValue = fnew;
            break;
            
        end
    end
    % Calculate x at K+1 iteration
    Xk      = Xk + alpha_K * DescentDirection_Xk;
    
    % THIS is Jacobian f ( x+alpha*dx)
    Jac_k  = [double(subs(Jac(1),{x1,x2},{Xk(1),Xk(2)}));
        double(subs(Jac(2),{x1,x2},{Xk(1),Xk(2)}))];
    
    % error definition
    err     = norm(Jac_k)/(1+abs(functionValue));
    Xold    = Xk;
    
    % update table
    cell_new       = {K,mat2str(Xold),mat2str(DescentDirection_Xk),...
        alpha_K,mat2str(Xk),err,functionValue};
    
    t_new = cell2table(cell_new);
    t_new.Properties.VariableNames = {'iteration','Initial Point', ...
        'Search Direction', 'Step length','new point','error',...
        'function value'};
    t = [t ; t_new];
    
    
    if err < tol
        break
    end
    disp("******")
end
disp('Optimal points of the function are: ')
disp(Xold)
disp('final optimal step size is: ')
disp(alpha_K)
disp('summary of values')
t

writetable(t,FileName)