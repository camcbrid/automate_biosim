function yf = runODE_montecarlo(A,sysname,q)
%run ODE for biomolecular system using reduced models including resource 
%sharing for many different parameters
%input the weighted, signed adjacency matrix
%inputs go in rows following adjacency matrix. This makes a non-square
%matrix with more rows than cols
%activation should be in [1,inf]; repression in [-1,0]
%edge order is done by cols first (order that A(:) returns)

%defaults
if nargin < 3
    q = 50;                   %number of trials to run
    if nargin < 2
        sysname = 'default';    %for output titles
        if nargin < 1
            %adjacency matrix
            A = [0 0;
                0 0;
                1 0];
        end
    end
end

uniformon = false;         %if want all nodes/edges to have same parameters
pfun = @(A) params_dist(A,uniformon);   %parameter function handle

%settings
tfinal = 1000;      %simulate from 0 to tfinal
umax = 10;          %sweep through inputs from 0 to umax with a trianglewave input
%init
n = min(size(A));   %number of nodes
yf = zeros(n,q);    %init output vector
[~,prng] = pfun(A); %#ok<ASGLU>

%run q simulations with structure of A and parameters drawn from the
%parameter distribution in params_dist
for ii = 1:q

    %create function handles: edit this function to change form of the dynamics
    [funs,p] = makefuns(A,pfun);
    
    %setup for ODE
    xmax = (p.RNAP*p.Ribo*p.k1*p.k2./(p.delta1*p.delta2)).*p.DNA./(p.Kp.*p.K2);
    x0 = xmax(:).*rand(n,1);
    f_cnst = @(t,x) dynamicsbio_static(t,x,funs,umax);
    opts = odeset('AbsTol',1e-6,'RelTol',1e-4,'Jacobian',...
        @(t,x) dynamicsbio_jac(t,x,funs,umax));
    
    %run ode
    sol = ode15s(f_cnst,[0 tfinal],x0,opts);
    
    %output
    yf(:,ii) = sol.y(:,end);
end

%condition output
yf(yf < 1e-4) = 0;
yf = real(yf);

%plot Monte Carlo results
fig = figure(3); clf;
%scatterplot of endpoints
subplot(211)
semilogy(real(yf'),'.')
ylim([1e-4,inf])
xlabel('simulation number')
ylabel('state concentrations')
title([strrep(sysname,'_',' '),', n = ',num2str(n)])

%plot histogram
subplot(212);
for jj = 1:n
    histogram(log10(yf(jj,:)),'BinMethod','Scott',...
        'Normalization','probability','DisplayStyle','stairs')
    hold on;
end
hold off;
ylabel('Frequency')
xlabel('log_{10}(steady state concentration)')

%save settings and output
savefile = [strrep(sysname,' ','_'),'_n',num2str(n),'_q',num2str(q)];
cd figs
savefig(fig,[savefile,'.fig']);
cd ..
cd output
save([savefile,'.mat'],'yf','prng','A','sysname','q');
cd ..

close(fig);
