function output = runODE_static(A, pfun, u, ploton)
%OUTPUT = runODE_static(A, PFUN, U, PLOTON)
%Simulate ODE with one constant input or no input to the system using
%reduced models including resource sharing. A is the weighted, signed
%adjacency matrix. External inputs go in the last row (n+1)st. If A is
%square, there is no external input when using random parameters.
%Activation edges > 0; repression edges < 0. Edge order is done by cols
%first (order that A(:) returns). PFUN is the parameter function handle or
%parameter struct. U is the concentration of the external input. PLOTON is
%a bool, which, if true, plots are displayed.

%defaults
if nargin < 4
    ploton = false;
    if nargin < 3
        %constant external input
        u = 0;
        if nargin < 2
            %parameter function handle or parameter struct
            pfun = @params_dist;
            if nargin < 1
                %adjacency matrix
                A = zeros(3);
            end
        end
    end
end

tic

%settings
addpath parameters utility
tfinal = 50;            %simulate from 0 to tfinal
n = min(size(A));       %number of nodes

%create function handles: edit this function to change form of the dynamics
[funs,p] = makefuns(A,pfun);

%setup for ODE
if all(isfield(p,{'RNAP','Ribo','k1','k2','delta1','delta2','DNA','Kp','K2'}))
    xmax = (p.RNAP*p.Ribo*p.k1*p.k2./(p.delta1*p.delta2)).*p.DNA./(p.Kp.*p.K2);
else; xmax = 5*ones(n,1);
end
x0 = 10.^(log10(xmax(:)).*rand(n,1));
f_cnst = @(t,x) dynamicsbio_static(t,x,funs,u,1,1);
opts = odeset('AbsTol',1e-6,'RelTol',1e-4,'Stats','off','Vectorized','on',...
    'Jacobian',@(t,x) dynamicsbio_jac(t,x,funs,u));

%run ode
sol = ode15s(f_cnst,[0 tfinal],x0,opts);
%outputs
tout = linspace(0,tfinal,700)';
yout = deval(sol,tout);
yf = yout(:,end);
toc

%Fourier decomposition
fftout = fft(yout(:,1:end)');
Fs = 1/(tout(2)-tout(1));
freq2 = Fs*(0:length(tout)/2)/length(tout);
[~,ind] = max(abs(fftout(2:end/2,:)));
if ind > 1; disp('oscillations'); end

%output data
output = struct;
output.t = tout;
output.y = yout;
output.yf = yf;
output.fftout = fftout;
output.p = p;

%plot
if ploton
    figure(2); clf;
    %proteins vs time
    subplot(211);
    h1 = semilogy(tout(3:end),yout(:,3:end));
    title('Time plot')
    xlabel('time [hrs]')
    ylabel('protein concentration [nM]')
    set(h1,'linewidth',1.5)
    %plot of fft decomposition
    subplot(212);
    h2 = loglog(freq2(1:end-1),abs(fftout(1:end/2,:)/length(tout)));
    xlabel('frequency [cycles/hr]')
    ylabel('fft magnitude [nM]')
    set(h2,'linewidth',1.5)
    drawnow
end

%mesh plot of det of Jacobian landscape
if n == 2 && ploton
    %settings
    ngrid = 20;     %number of grid points in each dimension
    %bounds for the mesh
    y1 = yout(1,:);
    y2 = yout(2,:);
    x1max = log10(max(y1))+.1;
    x2max = log10(max(y2))+.1;
    x1min = log10(min(y1(y1 > 0)))-.1;
    x2min = log10(min(y2(y2 > 0)))-.1;
    %init
    detdfdx = zeros(ngrid);
    detpath = zeros(length(yout),1);
    x1vec = logspace(x1min,x1max,ngrid);
    x2vec = logspace(x2min,x2max,ngrid);
    [X1,X2] = meshgrid(x1vec,x2vec);
    %loop through each grid point, evaluating the det of Jacobian at each
    for ii = 1:ngrid
        for jj = 1:ngrid
            y = [x1vec(jj);x2vec(ii)];
            [~,dfdz] = f_cnst(tfinal,y);
            detdfdx(ii,jj) = det(dfdz);
        end
    end
    %trajectory path on det Jacobian landscape
    for k = 1:length(yout)
        [~,dfdz2] = f_cnst(tfinal,yout(:,k));
        detpath(k) = det(dfdz2) + .01;
    end
    
    %plot mesh output
    figure(3); clf;
    surf(X1,X2,detdfdx,'edgecolor','none')
    set(gca,'xscale','log','yscale','log','zscale','linear')
    xlabel('x_1 concentration [nM]')
    ylabel('x_2 concentration [nM]')
    zlabel('det(df/dx) []')
    hold on
    plot3(yout(1,end),yout(2,end),detpath(end),'kx',...
        yout(1,1),yout(2,1),detpath(1),'ko',...
        yout(1,:),yout(2,:),detpath,'linewidth',3)
    xlim(10.^[x1min,x1max])
    ylim(10.^[x2min,x2max])
    alpha(0.5)
    view([-57.12 58.96])
end
