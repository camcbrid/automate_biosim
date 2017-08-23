function output = runfminconbio(A,pfun,u)
%run fmincon for arbitrary network using reduced models including resource
%sharing to find multiple equilibrium points. Input the weighted, signed 
%adjacency matrix, A. Inputs go in rows after adjacency matrix. External 
%inputs go in the last row (n+1)st. If no (n+1)st row, there is no external
%input if using random parameters; activation edges > 0; repression edges 
%< 0. Edge order is done by cols first (order that A(:) returns).

if nargin < 3
    %constant external input
    u = 90;
    if nargin < 2
        %parameter function handle or parameter struct
        pfun = @paramsRepCasc;
        if nargin < 1
            A = [0 -1;
                0 0;
                -1 0];
        end
    end
end

%settings
addpath parameters utility
ploton = true;             %if want to display a figure
n = min(size(A));           %number of nodes

tic
%create functions and normalize adjency matrix
[funs,p] = makefuns(A,pfun);  %edit this function to change dynamics

%setup fmincon problem
f_cnst = @(x) dynamicsbio_fmincon(0,x,funs,u,1,1);
if all(isfield(p,{'RNAP','Ribo','k1','k2','delta1','delta2','DNA','Kp','K2'}))
    xmax = (p.RNAP*p.Ribo*p.k1*p.k2./(p.delta1*p.delta2)).*p.DNA./(p.Kp.*p.K2);
else
    xmax = 5*ones(n,1);
end
x0 = 10.^(log10(xmax(:)).*rand(n,1));
lb = zeros(n,1);
ub = inf*ones(n,1);
opts = optimoptions('fmincon','Display','off',...%'Algorithm','sqp',...
    'SpecifyObjectiveGradient',true);

problem = createOptimProblem('fmincon','objective', ...
    f_cnst,'x0',x0,'lb',lb, 'ub',ub,'options',opts);

%run global search problem
% gs = GlobalSearch;
% [x,fval,eflag,stats,sols] = run(gs,problem)
ms = MultiStart;
[x,fval,eflag,stats,sols] = run(ms,problem,20);

%output data
output.x = x;
output.fval = fval;
output.flg = eflag;
output.stats = stats;
output.solutions = sols;
output.p = p;

numsols = length(output.solutions);
[X{1:numsols}] = deal(output.solutions.X);
[Fval{1:numsols}] = deal(output.solutions.Fval);

toc

y = cell2mat(X);
%Jacobian at each solution point
J = cell(0);
for ii = 1:numsols
    J{ii} = dynamicsbio_jac(0,y(:,ii),funs,u,1,1);
    disp(J{ii})
    disp(det(J{ii}))
end

%plot
if ploton
    figure(1); clf;
    
    %proteins vs time
    h1 = plot3(y(1,:),y(2,:),y(3,:),'o');
    set(gca,'xscale','log','yscale','log','zscale','log')
    title('State space')
    xlabel('x_1')
    ylabel('x_2')
    set(h1,'linewidth',1.5)
    
    drawnow
end
