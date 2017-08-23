function output = runODE_mulambda(A,pfun,u)
%run ODE for cascade using reduced models including resource sharing to
%find bifurcation. Input the weighted, signed adjacency matrix, A. Inputs
%go in rows after adjacency matrix. External inputs go in the last row
%(n+1)st. If no (n+1)st row, there is no external input if using random
%parameters; activation edges > 0; repression edges < 0. Edge order is done
%by cols first (order that A(:) returns).

if nargin < 3
    %constant external input
    u = .07;
    if nargin < 2
        %parameter function handle or parameter struct
        pfun = @(A) paramsActCasc(A);
        if nargin < 1
            A = [0 1;
                0 0;
                1 0];
        end
    end
end

%settings
addpath parameters utility
ploton = true;      %if want to display a figure
tfinal = 10000;     %simulate from 0 to tfinal
n = min(size(A));   %number of nodes

tic
%create function handles: edit this function to change form of the dynamics
[funs,p] = makefuns(A,pfun);

%setup for ODE
if all(isfield(p,{'RNAP','Ribo','k1','k2','delta1','delta2','DNA','Kp','K2'}))
    xmax = (p.RNAP*p.Ribo*p.k1*p.k2./(p.delta1*p.delta2)).*p.DNA./(p.Kp.*p.K2);
else
    xmax = 5*ones(n,1);
end
x0 = 10.^(log10(xmax(:)).*rand(n,1));
f_sweep  = @(t,x) dynamicsbio_mulambda(t,x,funs,u,tfinal);
opts = odeset('AbsTol',1e-9,'RelTol',1e-9,'Stats','off','Vectorized','on',...
    'Jacobian',@(t,x) dynamicsbio_jac(t,x,funs,max(u)));

%run ode
sol = ode15s(f_sweep,[0 tfinal],x0,opts);

%outputs
tout = linspace(0,tfinal,1200)';
yout = deval(sol,tout);
muvec = tout/tfinal;
lambdavec = muvec;
toc

%Jacobian
detf = zeros(length(tout),1);
if n < 10
    for jj = 1:length(tout)
        J = dynamicsbio_jac(tout(jj),yout(:,jj),funs,u,muvec(jj),lambdavec(jj));
        detf(jj) = det(J);
    end
end

%output data
output = struct;
output.y = yout;
output.t = tout;
output.mu = muvec;
output.lambda = lambdavec;
output.detf = detf;
output.p = p;

%plot
if ploton
    figure(6); clf;
    
    %protein concentration vs time
    subplot(221)
    h1 = semilogy(tout(1:end),yout(:,1:end));
    title('Time plot')
    xlabel('time')
    ylabel('protein concentration')
    set(h1,'linewidth',1.5)
    
    %protein concentration vs input
    subplot(223);
    h2 = semilogy(muvec(1:end),yout(:,1:end));
    title('Input response plot')
    xlabel('u')
    ylabel('protein concentration')
    set(h2,'linewidth',1.5)
    
    subplot(2,2,[2,4]);
    h3 = plot(muvec(3:end),detf(3:end));
    xlabel('u')
    ylabel('det(df/dx)')
    set(gca,'yaxislocation','right');
    set(h3,'linewidth',1.5)
    
    %input/response plot and input vs det(df/dx) plots
    figure(2); clf;
    h2 = semilogy(muvec(1:end),yout(:,1:end));
    title('Input response plot')
    xlabel('\mu')
    ylabel('protein concentration')
    set(h2,'linewidth',1.5)
    
    yyaxis right
    h3 = plot(muvec(3:end),(detf(3:end)),'--');
    ylabel('det(df/dx)')
    set(h3,'linewidth',1.5)
end

%mesh plot of det of Jacobian landscape
if n == 2 && ploton
    
    %settings
    ngrid = 20;     %number of grid points in each dimension
    f_jac = @(t,x) dynamicsbio_jac(t,x,funs,u);
    %bounds for the mesh
    y1 = yout(1,:);
    y2 = yout(2,:);
    x1max = log10(max(y1))+.1;
    x2max = log10(max(y2))+.1;
    x1min = log10(min(y1(y1 > 0)))-.1;
    x2min = log10(min(y2(y2 > 0)))-.1;
    %init
    Z = zeros(ngrid);
    Z2 = zeros(length(yout),1);
    x1vec = logspace(x1min,x1max,ngrid);
    x2vec = logspace(x2min,x2max,ngrid);
    [X1,X2] = meshgrid(x1vec,x2vec);
    tic
    %loop through each grid point, evaluating the det of Jacobian at each
    for ii = 1:ngrid
        for jj = 1:ngrid
            y = [x1vec(jj);x2vec(ii)];
            dfdz = f_jac(tfinal,y);
            Z(ii,jj) = det(dfdz);
        end
    end
    toc
    %trajectory path on det Jacobian landscape
    for k = 1:length(yout)
        dfdz2 = f_jac(tfinal,yout(:,k));
        Z2(k) = det(dfdz2) + .01;
    end
    toc
    
    %surface plot
    figure(7); clf;
    h = surf(X1,X2,Z);
    set(gca,'xscale','log','yscale','log','zscale','linear')
    shading flat
    set(h,'edgecolor','none')
    xlabel('x_1')
    ylabel('x_2')
    zlabel('det(df/dx)')
    hold on
    %plot trajectory on the surface
    plot3(yout(1,end),yout(2,end),Z2(end),'kx',...
        yout(1,1),yout(2,1),Z2(1),'k^',...
        yout(1,:),yout(2,:),Z2,'linewidth',3)
    xlim(10.^[x1min,x1max])
    ylim(10.^[x2min,x2max])
    alpha(0.5)
    view(2)
end

