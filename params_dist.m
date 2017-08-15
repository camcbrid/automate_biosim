function [p,prng] = params_dist(A,uniformon)
%assign parameters for production to each edge and degradation parameters
%for each protein. Takes in the weighted adacjency matrix, outputs
%parameter struct for functions h, g, a, and L to be used by makefuns.m
%using a specified distribution (normal, lognormal, uniform, etc)

%defaults
if nargin < 2
    uniformon = false;
    if nargin < 1
        error('not enough inputs')
    end
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
k = 1;                      %init node counter
edge_weights = A(A ~= 0);   %extract weighting of each edge

%need to pick for all duples of input combinations
if m > 1
    %find number of 
    twocombos = zeros(n,1);
    for jj = 1:n
        twocombos(jj) = nchoosek(indvec(jj),2);
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
if ~uniformon
    p.DNA = rng2out(DNArng,'norm',n);   %total DNA
    p.K0 = rng2out(K0rng,'logn',m);     %TF/DNA binding cnst
    p.K02 = rng2out(K0rng,'logn',q);    %TF/DNA binding cnst for two TFs at once
    p.K1 = rng2out(K1rng,'logn',m);     %RNAP/(TF/DNA) binding cnst
    p.Kp = rng2out(Kprng,'logn',n);     %RNAP/DNA leaky binding cnst
    p.K2 = rng2out(K2rng,'logn',n);     %mRNA/ribosome binding cnst
    p.s = rng2out(nrng,'norm',m);       %coopertivity
    p.K3 = rng2out(K3rng,'logn',n);     %protein/protease binding cnst
else
    p.DNA = rng2out(DNArng,'norm')*ones(1,n);   %total DNA
    p.K0 = rng2out(K0rng,'logn')*ones(1,m);     %TF/DNA binding cnst
    p.K02 = rng2out(K0rng,'logn')*ones(1,q);    %TF/DNA binding cnst for two TFs at once
    p.K1 = rng2out(K1rng,'logn')*ones(1,m);     %RNAP/(TF/DNA) binding cnst
    p.Kp = rng2out(Kprng,'logn')*ones(1,n);     %RNAP/DNA leaky binding cnst
    p.K2 = rng2out(K2rng,'logn')*ones(1,n);     %mRNA/ribosome binding cnst
    p.s = rng2out(nrng,'norm')*ones(1,m);       %coopertivity
    p.K3 = rng2out(K3rng,'logn')*ones(1,n);     %protein/protease binding cnst
end

%production parameters
hparams = struct;
%resource demand
hparams.J = p.DNA./p.Kp.*(1+p.k1*p.RNAP./(p.K2*p.delta1));
%leakiness
hparams.T = p.k2*p.k1*p.RNAP*p.Ribo*p.DNA./(p.Kp.*p.K2.*p.delta1);
%coopertivitity
hparams.n = p.s;
%loop through each edge
for ii = 1:m
    if edge_weights(ii) > 0
        %activation; max: T*a/b; min: T
        hparams.a(ii) = max([p.Kp(k)/p.K1(ii),p.K1(ii)/p.Kp(k)])*...
            abs(edge_weights(ii))/p.K0(ii); %activation strength
        hparams.b(ii) = 1./p.K0(ii);        %1/binding constant
    else
        %repression; max: T; min: T*a/b
        hparams.a(ii) = min([p.Kp(k)/p.K1(ii),p.K1(ii)/p.Kp(k)])/...
            (abs(edge_weights(ii))*p.K0(ii)); %activation strength
        hparams.b(ii) = 1./p.K0(ii);        %1/binding constant
    end
    %advance the node counter
    if any(ii == nodevec)
        %if the promoter has non-exclusive TF binding
        if indvec(k) > 1
            %allocate space for each node's c's and d's
            z = [0;cumsum(twocombos)];  %cumulative number of pairs of edges for each node
            rinds = nchoosek(nodevec(k)+1:nodevec(k+1),2);  %indicies of pairs of each edge
            for r = 1:(z(k+1) - z(k))
                s = z(k) + r;   %output counter
                %assign parameters for multiple protein binding to promoter
                d = 1./(p.K02(s)*p.K0(rinds(r,1))) + 1./(p.K02(s)*p.K0(rinds(r,2)));
                hparams.c(s) = max([p.Kp(k)/p.K1(ii),p.K1(ii)/p.Kp(k)])*d;
                hparams.d(s) = d;
            end
        end
        %advance the node counter
        k = k + 1;
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
p.delta = p.delta2;     %dilution rate constant
p.n = n;                %number of nodes
p.m = m;                %number of edges
p.h = hparams;          %production parameter struct
p.a = hparams;          %production resource sharing struct
p.g = gparams;          %degradation resource sharing struct

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


function out = rng2out(p_rng,distname,n)
%generates a random value within range. p_rng is a vector containing the
%upper and lower bounds for the range, n is the length of output, uniform
%is a binary option that, if true, generates an array of length n with
%identical elements

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
    out(out < 0) = 0;
elseif strcmp(distname,'logn')
    %lognormal distribution (input range is log10(range))
    mu = mean(log(10.^p_rng));
    sigma = (log(10^(max(p_rng) - min(p_rng))))/3;
    out = lognrnd(mu,sigma,[1,n]);
    out(out < 0) = 0;
elseif strcmp(distname,'unif')
    %uniform distribution
    out = min(p_rng) + (max(p_rng) - min(p_rng))*rand(1,n);
else
    mu = mean(p_rng);
    sigma = (max(p_rng) - min(p_rng))/3;
    out = random(distname,mu,sigma,[1,n]);
    out(out < 0) = 0;
end
