%run simulations of resource sensor for 3 different modulules (names 2,3,and 4)
%then compare their performance when placed together and quantify errors
%from resource sensor vs actual performance.

%system settings
u = logspace(0,6,15)';     %input (induction) vector
%system handles (add D for unconnected nodes)
addpath utility parameters
%adjency matrices
A1 = diag(ones(3-1,1),1);
A2 = diag(ones(4-1,1),1);
A = {A1,A2};
[A,n] = augment_input(A);
m = length(A);      %number of modules

%init
Acombo = cell(0);   %all two-combinations of adjacency matricies
pcombo = cell(0);   %parameters for modules together
funs = cell(0);     %function handles
Q = cell(0);        %estimated resource usage by each circuit
S = cell(0);        %estimated sensitivity of each circuit
S2 = cell(0);       %estimated sensitivity of each circuit
Fhat = cell(0);     %estimated change in output activation
w0 = cell(0);       %resource preturbation
wgF = cell(0);      %resource preturbation by output
x = cell(0);        %output states for module seperately
x2 = cell(0);       %output states for modules together
Qcalc = cell(0);    %calculated Q at steady state
beta = cell(0);     %tolerance S*Q
beta2 = cell(0);    %tolerance S2*Q

%% initialize parameters
p = params_dist(A);     %create parameter structs for each module
p = cell_params(p);     %set cellular properties to the same in each struct

%package parameters for running modules together
matches = nchoosek(1:m,2);
for i = 1:size(matches,1)
    pcombo{i} = combine_pstructs(p{matches(i,1)},p{matches(i,2)});
    Acombo{i} = combine_adjacency(A{matches(i,1)},A{matches(i,2)});
end

%% run 5 experiments
%measure resource usage
tic
for j = 1:m
    [x{j},Q{j},S{j},S2{j},Fhat{j},w0{j},wgF{j}] = RSexpms2(A{j},p{j},u,50,false);
    funs{j} = makefuns(A{j},p{j});
end
toc

%% run all pairs of modules
tic
%loop across each combination of pairs of modules
for r = 1:size(matches,1)
    %loop across levels of induction
    for q = 1:length(u)
        x2{r}(:,q) = runRSdynamics2(Acombo{r},pcombo{r},u(q),false);
    end
end
toc

%% actual Q and Fhat
for k = 1:length(A)
    for h = 1:length(u)
        Qcalc{k}(h,1) = 1./funs{k}.a(x{k}.z4(h,:)',u(h))-1;
    end
end

%% calculate sensitivities and tolerances for comparison

for i = 1:size(matches,1)
    beta{2*i-1} = S{matches(i,1)}.*Q{matches(i,2)};
    beta{2*i} = S{matches(i,2)}.*Q{matches(i,1)};
    beta2{2*i-1} = S2{matches(i,1)}.*Q{matches(i,2)};
    beta2{2*i} = S2{matches(i,2)}.*Q{matches(i,1)};
end

%% plotting settings
lwidth = 1.5;                           %line width
fntsze = 12;                            %axis font size
legfntsz = 10;                          %legend font size
pos = [145.4 110.2 787.6 530.4];        %position of figs 6,7,8

%% plot resource usage and sensitivities vs input, u
figure(14); clf;
subplot(131); hold off;
h31a = semilogx(u,[Q{:}]); hold on;
ax31 = gca; ax31.ColorOrderIndex = 1;
h31b = semilogx(u,[Qcalc{:}],'--');
ylabel('Q'); xlabel('Input, u [nM]');

subplot(132);
h32a = semilogx(u,[S{:}]); hold on;
ax32 = gca; ax32.ColorOrderIndex = 1;
h32b = semilogx(u,[S2{:}],'--');
ylabel('S'); xlabel('Input, u [nM]');

subplot(133);
h33a = semilogx(u,[Fhat{:}]);
ylabel('$\hat{F}$','Interpreter','Latex'); xlabel('Input, u [nM]');
legend(split(num2str(1:m),'  '),'Location','Best')
ax33 = gca;

set([ax31,ax32,ax33],'Fontsize',14);
set([h31a;h31b;h32a;h32b;h33a],'linewidth',lwidth)


%% module performance and prediction errors
figure(15); clf;
%actual vs predicted outputs
%loop through each module
matches_T = matches';
for d = 1:m
    subplot(2,m,d);
    
    betainds = find(matches_T(:) == d);
    %gives indicies in y for module d
    %index of combination of module in which the 'd'th module appears
    combotrial = sort(mod(find(matches == d) - 1, size(matches,1)) + 1);
    %whether the module appears first or second in each pair of modules
    modind = mod(find(matches' == d) - 1, 2) + 1;
    %lists index of output node for both modules in combinations where the
    %dth module appears
    nodeind = zeros(length(combotrial), 2);
    nodeind2 = zeros(length(combotrial), 1);
    y = zeros(length(u),length(combotrial));
    for k = 1:length(combotrial)
        nodeind(k,:) = cumsum([n{matches(combotrial(k),:)}]);
        nodeind2(k) = nodeind(k,modind(k));         %node index in combination trial
        y(:,k) = x2{combotrial(k)}(nodeind2(k),:)'; %output of correct module
    end
    h51a = loglog(u,x{d}.G4); hold on;  %module alone
    h51b = loglog(u, y);                %actual modules together
    ax51 = gca; ax51.ColorOrderIndex = 2;
    h51c = loglog(u, x{d}.G4./(1-[beta{betainds}]),'--');   %predicted modules together
    ylabel('Module output [nM]'); xlabel('Input, u [nM]');
    
    %error plots
    subplot(2,m,d+m); hold on;
    ax52 = gca; ax52.ColorOrderIndex = 2;
    h52 = semilogx(u,((x{d}.G4./(1-[beta{betainds}]))./y)-1);
    set(gca,'xscale','log')
    ylabel('Relative error'); xlabel('Input, u [nM]');
    
    set([ax51,ax52],'Fontsize',14);
    set([h51a;h51b;h51c;h52],'linewidth',lwidth);
end








return


%% plot simplified vs full models: steady state transfer curves
figure(1); clf;
subplot(131);
h11 = loglog(u,x{1}.z4(:,1:3)); hold on
ax = gca;
ax.ColorOrderIndex = 1;
loglog(u,[xA1,xA2,xA3],'--');
ax.ColorOrderIndex = 1;
loglog(u,[xA4,xA5,xA6],':'); hold off
xlabel('u')
title('A')
legend('x_1','x_2','x_3','Location','Best')
ylabel('Concentration [nM]')
set(gca,'Fontsize',14)

subplot(132);
ax = gca;
h12 = loglog(u,x{2}.z4(:,1:3)); hold on
ax.ColorOrderIndex = 1;
loglog(u,[xB1,xB2,xB3],'--');
ax.ColorOrderIndex = 1;
loglog(u,[xB4,xB5,xB6],':'); hold off
xlabel('u')
title('B')
legend('x_1','x_2','x_3','Location','Best')
ylabel('Concentration [nM]')
set(gca,'Fontsize',14)

subplot(133);
h13 = loglog(u,x{3}.z4(:,1:3)); hold on
ax = gca;
ax.ColorOrderIndex = 1;
loglog(u,[xC1,xC2,xC3],'--');
ax.ColorOrderIndex = 1;
loglog(u,[xC4,xC5,xC6],':'); hold off
xlabel('u')
title('C')
legend('x_1','x_2','x_3','Location','Best')
ylabel('Concentration [nM]')
set(gca,'Fontsize',14)

set([h11,h12,h13],'linewidth',lwidth)


%% Archive

% pG = struct();
% pG.k1 = 200;        % [/hr] rate of transcription (195 to 275)
% pG.k2 = 300;        % [/hr] TL rate (180 to 315)
% pG.delta1 = 12;     % [/hr] mRNA dilution rate
% pG.delta2 = .5;     % [/hr] protein dilution rate
% pG.RNAP = 500;      % [nM] total concentration of RNAP [nM] (2000-10000)
% pG.Ribo = 1000;     % [nM] total concentration of ribosomes [nM] (5800)
% pG.K1 = {560,560,560}; % [nM] binding RNAP to DNA promoter constant (150-560)
% %pG.Kp = {5e6,1e6,1e6}; % basal RNAP/promoter binding constant
% pG.K2 = {1e5,1e5,1e5}; % [nM] mRNA with ribosome binding constant (~10^4)
% pG.n = {1,1,1};
% %circuit A
% pA = struct(pG);
% pA.DNA = {150,200,60};
% pA.K0 = {5e3,1.6e3,5e3};
% pA.Kp = {5e6,1e6,1e6};
% %circuit B
% pB = struct(pG);
% pB.DNA = {150,200,60};
% pB.K0 = {5e3,1e7,5e3};
% pB.Kp = {5e6,1.6e3,1e6};
% %circuit C
% pC = struct(pG);
% pC.DNA = {150,200,60};
% pC.K0 = {5e7,1e7,5e3};
% pC.Kp = {5e3,1.6e3,1e6};
