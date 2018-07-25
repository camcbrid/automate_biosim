function [pRS, prng] = params_RS(A)
%[P2, PRNG] = params_RS(A)
%Assign parameters for production to each edge and degradation parameters
%for each protein. Takes in the weighted adacjency matrix, outputs
%parameter struct for functions h, g, a, and L to be used by makefuns.m
%using a specified distribution (normal, lognormal, uniform, etc). Binding
%constants use lognormal distributions and all other constants (catalytic
%rates, dilutions, total concentrations of resources, and coopertivity) use
%normal distributions.
%
%A is the weighted, signed adjacency matrix for the system.

%defaults
if nargin < 1
    error('not enough inputs')
end

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

%setup
B = (A ~= 0);               %logical adjacency matrix
nodevec = [0,cumsum(sum(B,1))]; %indicies of nodes in the edge counter
indvec = sum(B,1);          %number of incoming connections to node i
n = min(size(A));           %number of nodes
m = sum(B(A ~= 0));         %number of edges
edge_weights = A(A ~= 0);   %extract weighting of each edge
twocombos = zeros(n,1);     %init vector for picking duples of input combos to each node

%need to pick for all duples of input combinations
if m > 1
    %find number of pairs of edges incident to node jj
    for jj = 1:n
        if indvec(jj) > 1
            twocombos(jj) = nchoosek(indvec(jj),2);
        end
    end
    %number of additonal parameters required:
    q = sum(twocombos);
else
    q = 0;
end

%pick random parameters within each range
p = struct;                             %output parameter struct
p.k1 = rng2out(k1rng,'norm');           %TX rate
p.k2 = rng2out(k2rng,'norm');           %TL rate
p.k3 = rng2out(k3rng,'norm');           %degradation rate by protease
p.delta1 = rng2out(delta1rng,'norm');   %mRNA dilution/degradation
p.delta2 = rng2out(delta2rng,'norm');   %protein dilution
p.RNAP = rng2out(RNAPrng,'norm');       %total RNAP
p.Ribo = rng2out(Riborng,'norm');       %total ribosomes
p.Ptot = rng2out(Ptotrng,'norm');       %total protease
p.DNA = rng2out(DNArng,'norm',n+2);     %total DNA
p.K0 = rng2out(K0rng,'logn',m);         %TF/DNA binding cnst
p.K02 = rng2out(K0rng,'logn',q);        %TF/DNA binding cnst for two TFs at once
p.K1 = rng2out(K1rng,'logn',m);         %RNAP/(TF/DNA) binding cnst
p.Kp = rng2out(Kprng,'logn',n+2);       %RNAP/DNA leaky binding cnst
p.K2 = rng2out(K2rng,'logn',n+2);       %mRNA/ribosome binding cnst
p.s = rng2out(nrng,'norm',m);           %coopertivity
p.K3 = rng2out(K3rng,'logn',n+2);       %protein/protease binding cnst
p.m = m;                                %number of edges
p.n = n;                                %number of nodes

%split parameters and assign p.h, p.a, and p.g for each experiment
pRS = cell(0);
pRS{1} = assignHAG(exp_split_paras(p,1),nodevec,edge_weights,indvec,twocombos);
pRS{2} = assignHAG(exp_split_paras(p,2),nodevec,edge_weights,indvec,twocombos);
pRS{3} = assignHAG(exp_split_paras(p,3),nodevec,edge_weights,indvec,twocombos);
pRS{4} = assignHAG(exp_split_paras(p,4),nodevec,edge_weights,indvec,twocombos);
pRS{5} = assignHAG(exp_split_paras(p,5),nodevec,edge_weights,indvec,twocombos);

%record ranges
prng = struct;
prng.k1rng = k1rng;             prng.k2rng = k2rng;
prng.delta1rng = delta1rng;     prng.delta2rng = delta2rng;
prng.DNArng = DNArng;           prng.RNAPrng = RNAPrng;
prng.Riborng = Riborng;         prng.nrng = nrng;
prng.K0rng = K0rng;             prng.K1rng = K1rng;
prng.Kprng = Kprng;             prng.K2rng = K2rng;
prng.k3rng = k3rng;             prng.Ptotrng = Ptotrng;
prng.K3rng = K3rng;


function out = rng2out(p_rng, distname, n)
%generates a random value within range. p_rng is a vector containing the
%upper and lower bounds for the range, n is the length of output. DISTNAME
%may be 'unif', 'norm', 'logn', or any valid input distribution name to
%the function random.

%defaults to scalar and normal distribution
if nargin < 3
    n = 1;
    if nargin < 2
        distname = 'unif';
    end
end

%draw random samples from specified distribution with mu = mean(rng), and
%sigma = (max(rng) - min(rng))/3;
if strcmp(distname,'norm')
    %normal distribution
    mu = mean(p_rng);
    sigma = (max(p_rng) - min(p_rng))/3;
    out = normrnd(mu,sigma,[1,n]);
    out = abs(out);
elseif strcmp(distname,'logn')
    %lognormal distribution (input range is log10(range))
    mu = mean(log(10.^p_rng));
    sigma = (log(10^(max(p_rng) - min(p_rng))))/3;
    out = lognrnd(mu,sigma,[1,n]);
    out = abs(out);
elseif strcmp(distname,'unif')
    %uniform distribution
    out = min(p_rng) + (max(p_rng) - min(p_rng))*rand(1,n);
else
    mu = mean(p_rng);
    sigma = (max(p_rng) - min(p_rng))/3;
    out = random(distname,mu,sigma,[1,n]);
    out = abs(out);
end


function p2 = exp_split_paras(p, expnum)
%splits parameters contained in the struct P according to EXPNUM. For
%expnum == 1, 2 or 3,

%copy universal fields
p2.k1 = p.k1;           %TX rate
p2.k2 = p.k2;           %TL rate
p2.k3 = p.k3;           %degradation rate by protease
p2.delta1 = p.delta1;   %mRNA dilution/degradation
p2.delta2 = p.delta2;   %protein dilution
p2.RNAP = p.RNAP;       %total RNAP
p2.Ribo = p.Ribo;       %total ribosomes
p2.Ptot = p.Ptot;       %total protease
p2.m = p.m;             %number of edges
p2.n = p.n;             %number of nodes
switch expnum
    case 1
        %constitutive CFP
        p2.DNA = p.DNA(end);        %total DNA
        p2.Kp = p.Kp(end);          %RNAP/DNA leaky binding cnst
        p2.K2 = p.K2(end);          %mRNA/ribosome binding cnst
        p2.s = [];                  %coopertivity
        p2.K3 = p.K3(end);          %protein/protease binding cnst
    case 2
        %constitutive RFP
        p2.DNA = p.DNA(end-1);      %total DNA
        p2.Kp = p.Kp(end-1);        %RNAP/DNA leaky binding cnst
        p2.K2 = p.K2(end-1);        %mRNA/ribosome binding cnst
        p2.s = [];                  %coopertivity
        p2.K3 = p.K3(end-1);        %protein/protease binding cnst
    case 3
        %constitutive CFP and RFP
        p2.DNA = p.DNA(end-1:end);  %total DNA
        p2.Kp = p.Kp(end-1:end);    %RNAP/DNA leaky binding cnst
        p2.K2 = p.K2(end-1:end);    %mRNA/ribosome binding cnst
        p2.s = [];                  %coopertivity
        p2.K3 = p.K3(end-1:end);    %protein/protease binding cnst
    case 4
        %normal module
        p2.DNA = p.DNA(1:end-2);    %total DNA
        p2.K0 = p.K0;               %TF/DNA binding cnst
        p2.K02 = p.K02;             %TF/DNA binding cnst for two TFs at once
        p2.K1 = p.K1;               %RNAP/(TF/DNA) binding cnst
        p2.Kp = p.Kp(1:end-2);      %RNAP/DNA leaky binding cnst
        p2.K2 = p.K2(1:end-2);      %mRNA/ribosome binding cnst
        p2.s = p.s;                 %coopertivity
        p2.K3 = p.K3(1:end-2);      %protein/protease binding cnst
    case 5
        %module with constitutive RFP
        p2.DNA = p.DNA(1:end-1);    %total DNA
        p2.K0 = p.K0;               %TF/DNA binding cnst
        p2.K02 = p.K02;             %TF/DNA binding cnst for two TFs at once
        p2.K1 = p.K1;               %RNAP/(TF/DNA) binding cnst
        p2.Kp = p.Kp(1:end-1);      %RNAP/DNA leaky binding cnst
        p2.K2 = p.K2(1:end-1);      %mRNA/ribosome binding cnst
        p2.s = p.s;                 %coopertivity
        p2.K3 = p.K3(1:end-1);      %protein/protease binding cnst
    otherwise
        %output the input
        p2 = p;
end


function p2 = assignHAG(p, nodevec, edge_weights, indvec, twocombos)

k = 1;                             	%init node counter
%production parameters
hparams = struct;
%resource demand
hparams.J = p.DNA./p.Kp.*(1 + p.k1*p.RNAP./(p.K2*p.delta1));
%leakiness
hparams.T = p.k2*p.k1*p.RNAP*p.Ribo*p.DNA./(p.Kp.*p.K2.*p.delta1);
%coopertivitity
hparams.n = p.s;

%loop through each edge
for ii = 1:p.m
    if all(isfield(p,{'Kp','K1','K0'}))
        if edge_weights(ii) > 0
            %activation; max: T*a/b; min: T
            hparams.a(ii) = max([p.Kp(k)/p.K1(ii), p.K1(ii)/p.Kp(k)])*...
                (abs(edge_weights(ii))/p.K0(ii));   %activation strength
            hparams.b(ii) = 1./p.K0(ii);            %1/binding constant
        else
            %repression; max: T; min: T*a/b
            hparams.a(ii) = min([p.Kp(k)/p.K1(ii), p.K1(ii)/p.Kp(k)])/...
                (abs(edge_weights(ii))*p.K0(ii));   %activation strength
            hparams.b(ii) = 1./p.K0(ii);            %1/binding constant
        end
        %advance the node counter
        if any(ii == nodevec)
            %if the promoter has non-exclusive TF binding
            if indvec(k) > 1
                %allocate space for each node's c's and d's
                z = [0; cumsum(twocombos)];  %cumulative number of pairs of edges for each node
                rinds = nchoosek(nodevec(k)+1:nodevec(k+1),2);  %indicies of pairs of each edge
                for r = 1:(z(k+1) - z(k))
                    s = z(k) + r;   %output counter
                    %assign parameters for multiple protein binding to promoter
                    d = 1./(p.K02(s)*p.K0(rinds(r,1))) + ...
                        1./(p.K02(s)*p.K0(rinds(r,2)));
                    hparams.c(s) = max([p.Kp(k)/p.K1(ii),...
                        p.K1(ii)/p.Kp(k)])*d;
                    hparams.d(s) = d;
                end
            end
            %advance the node counter
            k = k + 1;
        end
    end
end

%shared degradation parameters degradation applies to all nodes (set K to
%inf if not degradated) could define index sets to use if had multiple
%proteases present
gparams = struct;
gparams.k = p.k3;       %catalytic rate of degradation
gparams.Ptot = p.Ptot;  %total amount of protease
gparams.K = p.K3;       %binding constant for each node

%combining into output parameter struct
p2.delta = p.delta2;     %dilution rate constant
p2.n = p.n;              %number of nodes
p2.m = p.m;              %number of edges
p2.h = hparams;          %production parameter struct
p2.a = hparams;          %production resource sharing struct
p2.g = gparams;          %degradation resource sharing struct

