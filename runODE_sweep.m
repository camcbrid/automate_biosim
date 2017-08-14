%simulate ODE using reduced models including resource sharing and sweep
%through a range of input from 0 to umax

%input the weighted, signed adjacency matrix
%inputs go in rows following adjacency matrix. This makes a non-square
%matrix with more rows than cols
%activation should be in [1,inf]; repression in [-1,0]
%edge order is done by cols first (order that A(:) returns)
%external inputs go in the last row (n+1)st. If no (n+1)st row, there is no external input
A = [0 1 1 1;
    1 0 1 1;
    1 1 0 1;
    1 1 1 0;
    1 0 0 0];

ploton = true;                          %if want to display a figure
uniformon = true;                       %if want all nodes/edges to have same parameters
pfun = @(A) params_dist(A,uniformon);   %parameter function handle

tic
%settings
tfinal = 10000;     %simulate from 0 to tfinal
umax = 100;         %sweep through inputs from 0 to umax with a trianglewave input
%setup
B = (A ~= 0);       %logical adjacency matrix
n = min(size(A));   %number of nodes

%create function handles: edit this function to change form of the dynamics
[funs,p] = makefuns(A,pfun);

%setup for ODE
xmax = (p.RNAP*p.Ribo*p.k1*p.k2./(p.delta1*p.delta2)).*p.DNA./(p.Kp.*p.K2);
x0 = xmax(:).*rand(n,1);
f_sweep  = @(t,x) dynamicsbio_sweep(t,x,funs,umax,tfinal,1,1);
opts = odeset('AbsTol',1e-6,'RelTol',1e-4,'Stats','on','Vectorized','on',...
    'Jacobian',@(t,x) dynamicsbio_jac(t,x,funs,umax));

%run ode
sol = ode15s(f_sweep,[0 tfinal],x0,opts);

%outputs
tout = linspace(0,tfinal,700)';
yout = deval(sol,tout);
uvec = umax*2*tout/tfinal;
uvec(tout>tfinal/2) = 2*umax*(1 - tout(end/2+1:end)/tfinal);
toc

%plot
if ploton
    figure(1); clf;
    
    %protein concentration vs time
    subplot(211)
    h1 = semilogy(tout(1:end),yout(:,1:end));
    title('Time plot')
    xlabel('time')
    ylabel('protein concentration')
    set(h1,'linewidth',1.5)
    
    %protein concentration vs input
    subplot(212);
    h2 = semilogy(uvec(1:end),yout(:,1:end));
    title('Input response plot')
    xlabel('u')
    ylabel('protein concentration')
    set(h2,'linewidth',1.5)
    
end

%mesh plot of det of Jacobian landscape
if n == 2 && ploton
    
    %settings
    ngrid = 15;     %number of grid points in each dimension
    f_cnst = @(t,x) dynamicsbio_static(t,x,funs,umax);
    %bounds for the mesh
    x1max = 1.2*log10(max(yout(1,:)));
    x2max = 1.2*log10(max(yout(2,:)));
    x1min = log10(min(y1(y1 > 0)));
    x2min = log10(min(y2(y2 > 0)));
    %init
    Z = zeros(ngrid);
    Z2 = zeros(length(yout),1);
    x1vec = logspace(x1min,x1max,ngrid);
    x2vec = logspace(x2min,x2max,ngrid);
    [X1,X2] = meshgrid(x1vec,x2vec);
    %loop through each grid point, evaluating the det of Jacobian at each
    for ii = 1:ngrid
        for jj = 1:ngrid
            y = [x1vec(jj);x2vec(ii)];
            [~,dfdz] = f_cnst(tfinal,y);
            Z(ii,jj) = log10(abs(det(dfdz)));
        end
    end
    %trajectory path on det Jacobian landscape
    for k = 1:length(yout)
        [~,dfdz2] = f_cnst(tfinal,yout(:,k));
        Z2(k) = log10(abs(det(dfdz2))) + .01;
    end
    
    %plot
    figure(3); clf;
    mesh(X1,X2,Z)
    set(gca,'xscale','log','yscale','log','zscale','linear')
    xlabel('x_1')
    ylabel('x_2')
    zlabel('log_{10}(det(df/dx))')
    hold on
    plot3(yout(1,end),yout(2,end),Z2(end),'kx',...
        yout(1,1),yout(2,1),Z2(1),'k^',...
        yout(1,:),yout(2,:),Z2,'linewidth',3)
    view(2)
end

