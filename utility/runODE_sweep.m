function output = runODE_sweep(A, pfun, urng, ploton)
%OUTPUT = runODE_sweep(A, PFUN, URNG, PLOTON)
%Run ODE sweeping external input, u, using reduced models including
%resource sharing. A is the weighted, signed adjacency matrix. External
%inputs go in the last row (n+1)st. If A is square, there is no external
%input when using random parameters. Activation edges > 0; repression edges
%< 0. Edge order is done by cols first (order that A(:) returns). PFUN is
%the parameter function handle or parameter struct. URNG is a vector with
%first element the minimum value of the desired input and second element as
%the maximum. PLOTON is a bool whether to display plots

%defaults
if nargin < 4
    ploton = true;
    if nargin < 3
        %sweep through inputs from 0 to umax with a trianglewave input
        urng = [86,94];
        if nargin < 2
            %parameter function handle or parameter struct
            pfun = @(A) paramsRepCasc(A);
            if nargin < 1
                A = [0 -1;
                    0 0;
                    -1 0];
            end
        end
    end
end

%settings
addpath parameters utility
tfinal = 10000;     %simulate from 0 to tfinal
n = min(size(A));   %number of nodes

tic
%create function handles: edit this function to change form of the dynamics
[funs,p] = makefuns(A,pfun);

%setup for ODE
if all(isfield(p,{'RNAP','Ribo','k1','k2','delta1','delta2','DNA','Kp','K2'}))
    xmax = (p.RNAP*p.Ribo*p.k1*p.k2./(p.delta1*p.delta2)).*p.DNA./(p.Kp.*p.K2);
else; xmax = 5*ones(n,1);
end
x0 = 10.^(log10(xmax(:)).*rand(n,1));
f_sweep = @(t,x) dynamicsbio_sweep(t,x,funs,urng,tfinal,1,1);
f_sweep2 = @(t,x) dynamicsbio_sweep(t,x,funs,urng,tfinal,0,0);
opts = odeset('AbsTol',1e-9,'RelTol',1e-9,'Stats','off','Vectorized','on',...
    'Jacobian',@(t,x) dynamicsbio_jac(t,x,funs,max(urng)));

%run ode
sol = ode15s(f_sweep,[0 tfinal],x0,opts);   %with resource sharing
sol2 = ode15s(f_sweep2,[0 tfinal],x0,opts); %no resource sharing
%outputs
tout = linspace(0,tfinal,1200)';
yout = deval(sol,tout);
yout2 = deval(sol2,tout);
uvec = (max(urng) - min(urng))*2*tout/tfinal + min(urng);
uvec(tout>tfinal/2) = 2*(max(urng) - min(urng))*...
    (1 - tout(end/2+1:end)/tfinal) + min(urng);
toc

%Jacobian
detf = zeros(length(tout),1);
if n < 10
    for jj = 1:length(tout)
        J = dynamicsbio_jac(tout(jj),yout(:,jj),funs,uvec(jj),1,1);
        detf(jj) = det(J);
    end
end

%output data
output = struct;
output.t = tout;
output.y = yout;
output.y2 = yout2;
output.u = uvec;
output.detf = detf;
output.p = p;

%plot
if ploton
    figure(4); clf;
    %protein concentration vs time
    subplot(221)
    h1 = semilogy(tout(3:end),yout(:,3:end));
    title('Time plot')
    xlabel('time [hrs]')
    ylabel('protein concentration [nM]')
    set(h1,'linewidth',1.5)
    %protein concentration vs input
    subplot(223);
    h2 = semilogy(uvec(3:end),yout(:,3:end));
    title('Input response plot')
    xlabel('u [nM]')
    ylabel('protein concentration [nM]')
    set(h2,'linewidth',1.5)
    %determinant of Jacobian vs input
    subplot(2,2,[2,4]);
    h3 = plot(uvec(3:end),detf(3:end));
    xlabel('u [nM]')
    ylabel('det(df/dx)')
    set(gca,'yaxislocation','right');
    set(h3,'linewidth',1.5)
    
    %input/response plot and input vs det(df/dx) plots
    figure(5); clf;
    h2 = semilogy(uvec(3:end),yout(:,3:end));
    hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    h3 = semilogy(uvec(3:end/2),yout2(:,3:end/2),'--');
    %     title('Input response plot')
    xlabel('u [nM]')
    ylabel('Concentration [nM]')
    set([h2,h3],'linewidth',1.5)
    %     h4 = semilogy(uvec(3:end),yout(1,3:end));
    %     hold on
    %     h5 = semilogy(uvec(3:end/2),yout2(1,3:end/2),'--');
    %     xlabel('u [nM]')
    %     ylabel('Concentration [nM]')
    %     set([h4,h5],'linewidth',1.5)
    %     yyaxis right
    %     h3 = plot(uvec(3:end),detf(3:end),'--');
    %     ylabel('det(df/dx)')
    %     set(h3,'linewidth',1.5)
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
            dfdz = dynamicsbio_jac(tfinal,y,funs,uvec(end),1,1);
            detdfdx(ii,jj) = det(dfdz);
        end
    end
    %trajectory path on det Jacobian landscape
    for k = 1:length(yout)
        dfdz2 = dynamicsbio_jac(tfinal,yout(:,k),funs,uvec(end),1,1);
        detpath(k) = det(dfdz2) + .01;
    end
    %plot
    figure(6); clf;
    surf(X1,X2,detdfdx,'edgecolor','none')
    set(gca,'xscale','log','yscale','log','zscale','linear')
    xlabel('x_1 concentration [nM]')
    ylabel('x_2 concentration [nM]')
    zlabel('det(df/dx) []')
    hold on
    plot3(yout(1,end),yout(2,end),detpath(end),'kx',...
        yout(1,3),yout(2,3),detpath(3),'ko',...
        yout(1,3:end),yout(2,3:end),detpath(3:end),'linewidth',3)
    xlim(10.^[x1min,x1max])
    ylim(10.^[x2min,x2max])
    alpha(0.5)
    view(2)
end