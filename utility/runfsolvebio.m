function output = runfsolvebio(A,pfun,u)
%run fsolve for arbitrary network using reduced models including resource
%sharing. Input the weighted, signed adjacency matrix, A. Inputs go in rows
%after adjacency matrix. External inputs go in the last row (n+1)st. If no 
%(n+1)st row, there is no external input if using random parameters; 
%activation edges > 0; repression edges < 0. Edge order is done by cols 
%first (order that A(:) returns).

if nargin < 3
    %constant external input
    u = 0;
    if nargin < 2
        %parameter function handle or parameter struct
        pfun = @(A) paramsARosc(A);
        if nargin < 1
            A = [1 1;
                -1 -1];
        end
    end
end

%settings
addpath parameters utility
n = min(size(A));       %number of nodes

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
f_cnst = @(x) dynamicsbio_static(0,x,funs,u);
fsolveopts = optimoptions(@fsolve,'Display','final','SpecifyObjectiveGradient',...
    true,'Algorithm','levenberg-marquardt');

%solve using fsolve
[x,fval,flg,soldata,jac] = fsolve(f_cnst,x0,fsolveopts);
if flg <= 0
    tfinal = 100;
    sol = ode15s(@(t,x) dynamicsbio_static(t,x,funs,u),[0,tfinal],x0);
    x2 = sol.y(:,end);
    [x,fval,flg,soldata,jac] = fsolve(f_cnst,x2,fsolveopts);
end
toc

%output data
output = struct;
output.x = x;
output.fval = fval;
output.flg = flg;
output.soldata = soldata;
output.jac = jac;
output.p = p;
