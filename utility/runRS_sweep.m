%run ODE sweeping external input, u, using reduced models including
%resource sharing. Input the weighted, signed adjacency matrix, A. Inputs go
%in rows after adjacency matrix. External inputs go in the last row
%(n+1)st. If no (n+1)st row, there is no external input if using random
%parameters; activation edges > 0; repression edges < 0. Edge order is done
%by cols first (order that A(:) returns).

%sweep through inputs from 0 to umax with a trianglewave input
urng = [0,50];
%parameter function handle or parameter struct
pfun1 = @(A) params_RS(A,1);
pfun2 = @(A) params_RS(A,2);
pfun3 = @(A) params_RS(A,3);
pfun4 = @(A) params_RS(A,4);
pfun5 = @(A) params_RS(A,5);
A = [0 1 0;
    0 0 1;
    0 0 0;
    1 0 0];

%settings
addpath parameters utility
ploton = true;      %if want to display a figure
tfinal = 1000;      %simulate from 0 to tfinal
n = min(size(A));   %number of nodes

A1 = [0];
A2 = [0];
A3 = [0 0;
    0 0];
A4 = A;
if max(size(A)) > n
    %external inputs
    A5 = [A(1:n,1:n), zeros(n,2);zeros(2,n+2);A(n+1:end,:),zeros(1,2)];
else
    %no external input
    A5 = [A, zeros(size(A,1),2);zeros(2,n+2)];
end

tic
%create function handles: edit this function to change form of the dynamics
[funs1,p1] = makefuns(A1,pfun);
[funs2,p2] = makefuns(A2,pfun);
[funs3,p3] = makefuns(A3,pfun);
[funs4,p4] = makefuns(A4,pfun);
[funs5,p5] = makefuns(A5,pfun);

%setup for ODE
if all(isfield(p5,{'RNAP','Ribo','k1','k2','delta1','delta2','DNA','Kp','K2'}))
    xmax = (p5.RNAP*p5.Ribo*p5.k1*p5.k2./(p5.delta1*p5.delta2)).*p5.DNA./(p5.Kp.*p5.K2);
else
    xmax = 5*ones(n,1);
end
x0 = 10.^(log10(xmax(:)).*rand(n+2,1));
x01 = x0(n+1);
x02 = x0(n+2);
x03 = x0(n+1:n+2);
x04 = x0(1:n);
x05 = x0(1:n+2);
f_sweep1 = @(t,x) dynamicsbio_sweep(t,x,funs1,urng,tfinal,1,1);
f_sweep2 = @(t,x) dynamicsbio_sweep(t,x,funs2,urng,tfinal,1,1);
f_sweep3 = @(t,x) dynamicsbio_sweep(t,x,funs3,urng,tfinal,1,1);
f_sweep4 = @(t,x) dynamicsbio_sweep(t,x,funs4,urng,tfinal,1,1);
f_sweep5 = @(t,x) dynamicsbio_sweep(t,x,funs5,urng,tfinal,1,1);
opts1 = odeset('AbsTol',1e-9,'RelTol',1e-9,'Stats','off','Vectorized','on',...
    'Jacobian',@(t,x) dynamicsbio_jac(t,x,funs1,max(urng)));
opts2 = odeset('AbsTol',1e-9,'RelTol',1e-9,'Stats','off','Vectorized','on',...
    'Jacobian',@(t,x) dynamicsbio_jac(t,x,funs2,max(urng)));
opts3 = odeset('AbsTol',1e-9,'RelTol',1e-9,'Stats','off','Vectorized','on',...
    'Jacobian',@(t,x) dynamicsbio_jac(t,x,funs3,max(urng)));
opts4 = odeset('AbsTol',1e-9,'RelTol',1e-9,'Stats','off','Vectorized','on',...
    'Jacobian',@(t,x) dynamicsbio_jac(t,x,funs4,max(urng)));
opts5 = odeset('AbsTol',1e-9,'RelTol',1e-9,'Stats','off','Vectorized','on',...
    'Jacobian',@(t,x) dynamicsbio_jac(t,x,funs5,max(urng)));

%run ode
sol1 = ode15s(f_sweep1,[0 tfinal],x01,opts1);
sol2 = ode15s(f_sweep2,[0 tfinal],x02,opts2);
sol3 = ode15s(f_sweep3,[0 tfinal],x03,opts3);
sol4 = ode15s(f_sweep4,[0 tfinal],x04,opts4);
sol5 = ode15s(f_sweep5,[0 tfinal],x05,opts5);

%outputs
tout = linspace(0,tfinal,1200)';
yout1 = deval(sol1,tout);
yout2 = deval(sol2,tout);
yout3 = deval(sol3,tout);
yout4 = deval(sol4,tout);
yout5 = deval(sol5,tout);
uvec = (max(urng) - min(urng))*2*tout/tfinal + min(urng);
uvec(tout>tfinal/2) = 2*(max(urng) - min(urng))*...
    (1 - tout(end/2+1:end)/tfinal) + min(urng);
toc

%Jacobian
detf = zeros(length(tout),1);
if n < 10
    for jj = 1:length(tout)
        J = dynamicsbio_jac(tout(jj),yout4(:,jj),funs4,uvec(jj),1,1);
        detf(jj) = det(J);
    end
end

%output data
output = struct;
output.t = tout;
output.y = yout4;
output.u = uvec;
output.detf = detf;
output.p = p5;

%plot
if ploton
    figure(4); clf;
    %protein concentration vs time
    subplot(221)
    h1 = semilogy(tout(3:end),yout4(:,3:end));
    title('Time plot')
    xlabel('time [hrs]')
    ylabel('protein concentration [nM]')
    set(h1,'linewidth',1.5)
    
    %protein concentration vs input
    subplot(223);
    h2 = semilogy(uvec(3:end),yout4(:,3:end));
    title('Input response plot')
    xlabel('u [nM]')
    ylabel('protein concentration [nM]')
    set(h2,'linewidth',1.5)
    %input vs det(df/dx) plot
    subplot(2,2,[2,4]);
    h3 = plot(uvec(3:end),detf(3:end));
    xlabel('u [nM]')
    ylabel('det(df/dx)')
    set(gca,'yaxislocation','right');
    set(h3,'linewidth',1.5)
    
    %input/response plot and input vs det(df/dx) plot
    figure(5); clf;
    h2 = semilogy(uvec(3:end),yout4(:,3:end));
    title('Input response plot')
    xlabel('u [nM]')
    ylabel('protein concentration [nM]')
    set(h2,'linewidth',1.5)
    yyaxis right
    h3 = plot(uvec(3:end),detf(3:end),'--');
    ylabel('det(df/dx)')
    set(h3,'linewidth',1.5)
end


%mesh plot of det of Jacobian landscape
if n == 2 && ploton
    %settings
    ngrid = 20;     %number of grid points in each dimension
    %bounds for the mesh
    y1 = yout4(1,:);
    y2 = yout4(2,:);
    x1max = log10(max(y1))+.1;
    x2max = log10(max(y2))+.1;
    x1min = log10(min(y1(y1 > 0)))-.1;
    x2min = log10(min(y2(y2 > 0)))-.1;
    
    %init
    detdfdx = zeros(ngrid);
    detpath = zeros(length(yout4),1);
    x1vec = logspace(x1min,x1max,ngrid);
    x2vec = logspace(x2min,x2max,ngrid);
    [X1,X2] = meshgrid(x1vec,x2vec);
    %loop through each grid point, evaluating the det of Jacobian at each
    for ii = 1:ngrid
        for jj = 1:ngrid
            y = [x1vec(jj);x2vec(ii)];
            dfdz = dynamicsbio_jac(tfinal,y,funs4,uvec(end),1,1);
            detdfdx(ii,jj) = det(dfdz);
        end
    end
    %trajectory path on det Jacobian landscape
    for k = 1:length(yout4)
        dfdz2 = dynamicsbio_jac(tfinal,yout4(:,k),funs4,uvec(end),1,1);
        detpath(k) = det(dfdz2) + .01;
    end
    
    %plot
    figure(12); clf;
    surf(X1,X2,detdfdx,'edgecolor','none')
    set(gca,'xscale','log','yscale','log','zscale','linear')
    xlabel('x_1 concentration [nM]')
    ylabel('x_2 concentration [nM]')
    zlabel('det(df/dx) []')
    hold on
    plot3(yout4(1,end),yout4(2,end),detpath(end),'kx',...
        yout4(1,3),yout4(2,3),detpath(3),'ko',...
        yout4(1,3:end),yout4(2,3:end),detpath(3:end),'linewidth',3)
    xlim(10.^[x1min,x1max])
    ylim(10.^[x2min,x2max])
    alpha(0.5)
    view(2)
end
