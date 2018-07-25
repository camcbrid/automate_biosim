function [y_out, figout, sol] = runRSdynamics2(A, p, u, ploton, fignum, linesty)
%[Y_OUT, FIGOUT, SOL] = runRSdynamics2(A, P, U, PLOTON, FIGNUM, LINESTY)
%run ode for resource_sensor_full2,3,4 (all 2-node activation/repression
%cascades) and output the steady states of all proteins in system and free
%RNAP/ribosomes.
%
%A is the adjacency matrix, P is the parameter struct returned by a
%parameter function, U is the concentration of an external input, PLOTON is
%a bool which, if true, plots figures, FIGNUM is the number of the desired
%figure to output the plot to, LINESTY is a string for the line style on
%the plot output. Y_OUT is the concentration of all states at steady state,
%FIGOUT is the figure handle of the plot output, SOL is the solution struct
%of the ODE simulation.

%defaults
if nargin < 4
    ploton = false;
    if nargin < 3
        u = 0;
        if nargin < 2
            if nargin < 1
                A = zeros(3);
            end
            p = params_dist(A);
        end
    end
end

n = size(A,2);
%create function handles for dynamics
funs = makefuns(A,p);

%init ode solver
tfinal = 200;
odefun = @(t,x) dynamicsbio_static(t,x,funs,u,1,1);
x0 = ones(n,1);
opts = odeset('AbsTol',1e-6,'RelTol',1e-5,'Stats','off');%,'Vectorized','on',...
%'Jacobian',@(t,x) dynamicsbio_jac(t,x,funs,u));

%run ode
sol = ode15s(odefun,[0 tfinal],x0,opts);
tout = linspace(0,tfinal,200)';
yout = deval(sol,tout)';
y_out = yout(end,:);

%plot output
if ploton
    if nargin >= 5
        fig = figure(fignum); hold on;
        ax = gca; ax.ColorOrderIndex = 1;
        if nargin < 6
            linesty = '-';
        end
    else
        fig = figure; clf;
        linesty = '-';
    end
    %proteins
    h1 = semilogy(tout,yout,linesty);
    ylabel('Protein Concentration [nM]')
    xlabel('Time [hrs]')
    legend(split(num2str(1:n),'  '),'Location','Best')
    set(h1,'linewidth',1.5)
end

%output figure handle if desired
if nargout > 1
    figout = gcf;
end