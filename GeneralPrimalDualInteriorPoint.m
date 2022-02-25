function GeneralPrimalDualInteriorPoint(A, b,c, Name)


%% varaiable definition
tol = 1e-6;
n = size(A,2) ;   % number of columns of A
m = size(A,1);   % number of rows of A
itr = 10;
etha = 0.95;

% set output directories
path = fullfile('.','tables');

if ~exist(path, 'dir')
    mkdir(path)
end

%% function definition
alpha = @(v, dv) min (1, min(-v(dv<0)./dv(dv<0)));
y = @(x, z, n) x'*z/n;
yaff = @(x, ax, dx, z, az, dz, n )  (x + ax* dx)'*(z+az*dz)/n;
tau = @(yaff, y) (yaff/y)^3;
rp = @(A,x,b) -(A*x -b);
rd = @(A,pi,z,c) -(A' * pi + z -c);


set_d_aff = @(A,rp,rd,x,z,n,m)  ...
    [  zeros(n) A' eye(n);
    A zeros(m) zeros(m,n);
    diag(z) zeros(n,m) diag(x)] \[rd ; rp ; diag(-diag(x)*diag(z))];

set_d = @(A,rp,rd,x,z,n,m, dx_aff, dz_aff, tk, yk)  ...
    [    zeros(n) A' eye(n);
    A zeros(m) zeros(m,n);
    diag(z) zeros(n,m) diag(x)]\[rd ; rp ;
    diag(-diag(x)*diag(z)-diag(dx_aff)*diag(dz_aff)+tk*yk*eye(n))];


dx = @(d, n) d(1:n);
dpi = @(d,n,m) d(n+1:n+m);
dz = @(d,n,m) d(n+m+1:end);


stopConditions = @(A,x,b,pi,z,c) [norm(A*x-b);
                                            norm(A'*pi+z-c);
                                            x'*z];


%% initial point generator
c
[x0,z0,pi0] = initialPointGeneration(A, b, c);
sc0 = stopConditions(A,x0,b,pi0,z0,c);
t = table(0,x0(1),x0(2),x0(3),x0(4), c'*x0, pi0(1),pi0(2),...
    z0(1),z0(2),z0(3),z0(4), b'*pi0, 0 ,sc0(1), sc0(2), sc0(3));
t.Properties.VariableNames = ...
    {'K','x1','x2','x3','x4','cTx','\pi1','\pi2','z1','z2','z3','z4','bTpi',...
    '\tau','||Ax-b||','||A\pi+z-c||','xz'};

%% Iteration generation

for K = 1:itr
   
    % CALCULATING RESIDUALS
    rp0 = rp(A,x0,b);
    rd0 = rd(A,pi0,z0,c);
    
    
    % STEP 1
    d_aff = set_d_aff(A,rp0,rd0,x0,z0,n,m);
    dx_aff_0 = dx(d_aff, n);
    dpi_aff_0 = dpi (d_aff,n,m);
    dz_aff_0 = dz (d_aff,n,m);
    
    alpha_aff_0_x = alpha(x0,dx_aff_0);
    alpha_aff_0_z = alpha(z0,dz_aff_0);
    
    y0 = y(x0,z0,n);
    
    y0_aff = yaff (x0,alpha_aff_0_x,dx_aff_0,z0,alpha_aff_0_z,dz_aff_0,n);
    
    t0 = tau(y0_aff,y0);
    
    
    d = set_d(A,rp0,rd0,x0,z0,n,m, dx_aff_0, dz_aff_0, t0, y0);
    
    
    dx_0 = dx(d, n);
    dpi_0 = dpi (d,n,m);
    dz_0 = dz (d,n,m);
    
    
    alpha_max_x_0 = alpha(x0,dx_0);
    alpha_max_z_0 = alpha(z0,dz_0);
    
    alpha_x_0 = min (1, etha*alpha_max_x_0);
    alpha_z_0 = min (1, etha*alpha_max_z_0);
    
    
    xk = x0 + alpha_x_0 * dx_0;
    pik = pi0 + alpha_z_0 * dpi_0;
    zk = z0 + alpha_z_0 * dz_0;
    
    sc = stopConditions(A,xk,b,pik,zk,c);
    
    if all(sc<= tol)
        fprintf("tolerence has reached at iteration %d \n", K)
        break;
    end
    
    x0 = xk;
    pi0 = pik;
    z0 = zk;
    
    tk = table(K, xk(1),xk(2),xk(3),xk(4), ...
        c' * xk, pik(1),pik(2), zk(1),zk(2),zk(3),zk(4),...
        b'*pik, t0, sc(1),sc(2), sc(3));
    tk.Properties.VariableNames = ...
        {'K','x1','x2','x3','x4','cTx','\pi1','\pi2','z1','z2','z3','z4','bTpi',...
        '\tau','||Ax-b||','||A\pi+z-c||','xz'};
    t = [t;tk];
    
    if K == 6
        x0
        xk
        alpha_x_0 
        dx_0
    end
end


disp("saving results in: ");
FileName = fullfile(path,Name);
disp(FileName)
writetable(t ,FileName);



end





