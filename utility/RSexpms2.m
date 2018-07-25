function [x, Q, S, S2, Fhat, w0, wgF] = RSexpms2(A, p, u, RFPtot, ploton)
%[X, Q, S, S2, FHAT, W0, WGF] = RSexpms2(A, P, U, RFPTOT, PLOTON)
%Run five experiments to characterize resource usage and sensitivity of
%module With adjacency matrix A. Finds Q, S, and Fhat
%
%A is the adjacency matrix, P is the parameter struct, U is the
%concentration of the external input, RFPTOT is the concentration of the
%DNA of the distrubance node, PLOTON is a bool, which, if true, outputs
%plots. X is a struct containing all data for outputing, Q is the metric
%for resource usage by the system, S is the metric for sensitivity of the
%output of the module to resource disturbances using the left endpoint 
%method, S2 is the resource sensivity metric calculated the using the right
%endpoint method, FHAT is the ratio of activation of the output node with
%and without the RFP distrubance, W0 is the resources used by the RFP
%disturbance, and WGF is the resources used by the constitutive GFP.

%defaults
if nargin < 5
    ploton = false;
    if nargin < 4
        RFPtot = 30;
        if nargin < 3
            u = logspace(0,10,15);
            if nargin < 2
                if nargin < 1
                    A = ones(3) - eye(3);
                    A = [A;1,zeros(1,size(A,1)-1)];     %augment with input
                end
                p = params_dist(A);
            end
        end
    end
end

if iscell(p) && iscell(A)
    x = cell(0);
    for j = 1:length(A)
        x{j} = RSexpms2(A,p,u,RFPtot,ploton);
    end
    return
end

%initialize
n = min(size(A));
z4 = zeros(length(u),n);
z5 = zeros(length(u),n+1);

%create adjacency matricies for each experiment
A1 = zeros(1);                      %YFP alone
A2 = zeros(1);                      %RFP alone
A3 = zeros(2);                      %YFP and RFP together
A4 = A;                             %module alone with input
A5 = combine_adjacency(A4,zeros(1));%module with RFP

%parameters
p1 = params_dist(A1);               %YFP alone
p2 = params_dist(A2);               %RFP alone
p3 = combine_pstructs(p1,p2);       %YFP and RFP together
p4 = p;                             %module alone
p5 = combine_pstructs(p4,p2);       %module with RFP
p0 = cell_params(p4,p1,p2,p3,p5);   %set cellular properties to be the same
[p4,p1,p2,p3,p5] = deal(p0{:});
%set DNA copy numbers
p1.DNA = RFPtot/10;
p2.DNA = RFPtot;
p3.DNA = [RFPtot/10,RFPtot];
p5.DNA(end) = RFPtot;

sol4 = cell(0);
sol5 = cell(0);

%experiment 1: only GFP
[z1,fig1,sol1] = runRSdynamics2(A1,p1,0,ploton);
%experiment 2: only RFP
[z2,~,sol2] = runRSdynamics2(A2,p2,0,ploton,fig1,'--');
%experiment 3: only GFP and RFP
[z3,~,sol3] = runRSdynamics2(A3,p3,0,ploton,fig1,'-.');
for ii = 1:length(u)    %sweep through induction levels
    %experiment 4: module with reporter
    [z4(ii,:),fig2,sol4{ii}] = runRSdynamics2(A4,p4,u(ii),ploton);
    %experiment 5: module with reporter and RFP sensor
    [z5(ii,:),~,sol5{ii}] = runRSdynamics2(A5,p5,u(ii),ploton,fig2,'--');
end

%calculate desired quantities from 5 experimental results at steady state
Y1 = z1(end,1);                     R2 = z2(end,1);
Y3 = z3(end,1);                     R3 = z3(end,2);
G4 = z4(:,n);
G5 = z5(:,n);                       R5 = z5(:,n+1);
Ghat0 = Y3/Y1;                      Rhat0 = R3/R2;
Ghat = G5./G4;                      Rhat = R5/R2;

%resources used by RFP (scalar)
w0 = (1-Ghat0)/(Rhat0 + Ghat0 - 1);
%resources used by GFP (scalar)
wgF = (1-Rhat0)/(Rhat0 + Ghat0 - 1);
%change in induction of GFP (vector)
Fhat = (1+w0)*Ghat./(1+w0*(1-Rhat));

%resource measures (quantity and sensitivity)
Q = (1./Rhat - 1)*(1+w0);     %quantity of resources used by module
S = (Ghat-1)./(w0.*Fhat);     %sensitivity of output to resource disturbance
S2 = (Ghat-1)./(w0*Ghat.*Fhat);%sensitivity of output to resource disturbance

%put everything into a struct for output
x.G1 = Y1;              x.R2 = R2;
x.G3 = Y3;              x.R3 = R3;
x.G4 = G4;              x.G5 = G5;
x.R5 = R5;
x.Ghat0 = Ghat0;        x.Rhat0 = Rhat0;
x.Ghat = Ghat;          x.Rhat = Rhat;
x.z1 = z1;              x.z2 = z2;
x.z3 = z3;              x.z4 = z4;
x.z5 = z5;
x.w0 = w0;              x.wgF = wgF;
x.Fhat = Fhat;          x.Q = Q;
x.S = S;                x.S2 = S2;
x.A1 = A1;              x.A2 = A2;
x.A3 = A3;              x.A4 = A4;
x.A5 = A5;
x.p1 = p1;              x.p2 = p2;
x.p3 = p3;              x.p4 = p4;
x.p5 = p5;              x.p0 = p0;
x.sol1 = sol1;          x.sol2 = sol2;
x.sol3 = sol3;          x.sol4 = sol4;
x.sol5 = sol5;

%plot outputs
if ploton
    figure; clf;
    
    subplot(221);
    h1 = loglog(u,[Y1*ones(size(u)),Y3*ones(size(u)),R2*ones(size(u)),...
        R3*ones(size(u)),R5]);
    ylabel('YFP, GFP, RFP')
    legend('Y_1','Y_3','R_2','R_3','R_5','Location','Best')
    title('Experiments')
    set(gca,'Fontsize',14)
    
    subplot(222);
    h2 = loglog(u,[G4,G5,R5]);
    ylabel('YFP, GFP, RFP')
    legend('G_4','G_5','R_5','Location','Best')
    set(gca,'Fontsize',14)
    
    subplot(223);
    h3 = semilogx(u,[Ghat0*ones(size(u)),Rhat0*ones(size(u)),Ghat,Rhat]);
    ylabel('GR''')
    legend('G_0','R_0','G''','R''','Location','Best')
    set(gca,'Fontsize',14)
    
    subplot(224);
    h4 = semilogx(u,[Q,S]);
    ylabel('Q/S')
    legend('Q','S','Location','Best')
    set(gca,'Fontsize',14)
    
    set([h1;h2;h3;h4],'linewidth',1.5)
end
