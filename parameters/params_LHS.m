function p = params_LHS(A)
%assign parameters for production to each edge and degradation parameters
%for each protein. Takes in the weighted adacjency matrix, outputs
%parameter struct for functions h, g, a, and L to be used by makefuns.m

if nargin < 1
    error('not enough inputs')
end
B = (A ~= 0);                       %logical adjacency matrix
n = min(size(A));                   %number of nodes
m = sum(B(A ~= 0));                 %number of edges

%extract weighting of each edge
edge_weights = A(A ~= 0);

%output parameter struct
p = struct;

%set the ranges for each parameter
k1rng = [195,275];      % [/hr] rate of transcription (195 to 275)
k2rng = [180,315];      % [/hr] TL rate (180 to 315)
delta1rng = [7,20];     % [/hr] mRNA dilution rate (7 to 20)
delta2rng = [.4,1];     % [/hr] protein dilution rate (0.4 to 1)
DNArng = [5,200];       % [nM] copy numbers of DNA promoters (5 to 200)
RNAPrng = [500,2000];   % [nM] total concentration of RNAP (2000 to 10000)
Riborng = [100,1000];   % [nM] total concentration of ribosomes (5800)
nrng = [1,3];           % [] coopertivity (1-3)
%powers of 10
K0rng = [1,4];          % [nM] binding TF with DNA promoter constant ()??
K1rng = [3,6];          % [nM] binding RNAP to DNA promoter constant (150-2000)
Kprng = [3,6];          % [nM] basal RNAP/promoter binding constant (150-560)
K2rng = [3,6];          % [nM] mRNA with ribosome binding constant (~10^4)
%protease parameters
k3rng = [270,270];      % [/hr] catalytic rate of degradation by protease
Ptotrng = [0,0];        % [nM] total amount of protease
K3rng = log10([1200,1200]);    % [nM] protein with protease binding constant

%pick random parameters within each range
p.k1 = rng2out(k1rng);          %TX rate
p.k2 = rng2out(k2rng);          %TL rate
p.k3 = rng2out(k3rng);          %degradation rate by protease
p.delta1 = rng2out(delta1rng);  %mRNA dilution/degradation
p.delta2 = rng2out(delta2rng);  %protein dilution
p.RNAP = rng2out(RNAPrng);      %total RNAP
p.Ribo = rng2out(Riborng);      %total ribosomes
p.Ptot = rng2out(Ptotrng);      %total protease
p.DNA = rng2out(DNArng);        %total DNA
p.K0 = 10.^rng2out(K0rng);      %TF/DNA binding cnst
p.K1 = 10.^rng2out(K1rng);      %RNAP/(TF/DNA) binding cnst
p.Kp = 10.^rng2out(Kprng);      %RNAP/DNA leaky binding cnst
p.K2 = 10.^rng2out(K2rng);      %mRNA/ribosome binding cnst
p.s = rng2out(nrng);            %coopertivity
p.K3 = rng2out(K3rng);          %protein/protease binding cnst

%production parameters
hparams = struct;
hparams.J2 = num2cell(0*ones(1,n));
%resource demand
hparams.J = num2cell(p.DNA./p.Kp.*(1+p.k1*p.RNAP./(p.K2*p.delta1)));
%leakiness
hparams.T = num2cell(p.k2*p.k1*p.RNAP*p.Ribo*p.DNA./(p.Kp.*p.K2.*p.delta1));
%coopertivitity
hparams.n = num2cell(p.s);

%loop through each edge
for ii = 1:m
    if edge_weights(ii) > 0
        %activation; max: T*a/b; min: T
        %activation strength
        hparams.a{ii} = (max([p.Kp/p.K1,p.K1/p.Kp])/p.K0)*abs(edge_weights(ii));
        hparams.b{ii} = 1./p.K0;                       %1/binding constant
        %disp('activation')
    else
        %repression; max: T; min: T*a/b
        %activation strength
        hparams.a{ii} = (min([p.Kp/p.K1,p.K1/p.Kp])/p.K0)*abs(edge_weights(ii));
        hparams.b{ii} = 1./p.K0;                       %1/binding constant
        %disp('repression')
    end
end

%shared degradation parameters degradation applies to all nodes (set K to
%inf if not degradated) could define index sets to use if had multiple
%proteases present
gparams = struct;
gparams.k = num2cell(p.k3);         %catalytic rate of degradation
gparams.Ptot = num2cell(p.Ptot);    %total amount of protease
gparams.K = num2cell(p.K3);         %binding constant for each node

%combining into output parameter struct
p.delta = p.delta2;             %dilution rate constant
p.n = n;                        %number of nodes
p.m = m;                        %number of edges
p.h = hparams;                  %production parameter struct
p.a = hparams;                  %production resource sharing struct
p.g = gparams;                  %degradation resource sharing struct


function out = rng2out(p_rng,n,distname)
%generates a random value within range. p_rng is a vector containing the
%upper and lower bounds for the range, n is the length of output, uniform
%is a binary option that, if true, generates an array of length n with
%identical elements

%defaults
if nargin < 3
    distname = 'logn';
    if nargin < 2
        n = 1;
    end
end

if strcmp(distname,'norm')
    mu = mean(p_rng);
    sigma = (max(p_rng) - min(p_rng))/2;
    out = normrnd(mu,sigma,[n,1]);
elseif ~strcmp(distname,'norm')
    mu = mean(p_rng);
    sigma = (max(p_rng) - min(p_rng))/2;
    out = random(distname,mu,sigma,[n,1]);
end
