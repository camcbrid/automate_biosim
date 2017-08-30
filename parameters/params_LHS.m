function [p,prng] = params_LHS(A,numsamples)
%assign parameters for production to each edge and degradation parameters
%for each protein. Takes in the weighted adacjency matrix, outputs
%parameter struct for functions h, g, a, and L to be used by makefuns.m
%using a specified distribution (normal, lognormal, uniform, etc)

%defaults
if nargin < 2
    numsamples = 100;
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
edge_weights = A(A ~= 0);   %extract weighting of each edge

%need to pick for all duples of input combinations
if m > 1
    %find number of pairs of edges incident to node jj
    twocombos = zeros(n,1);
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

%create latin hypercube
X = lhsdesign(numsamples, 8 + 4*n + 3*m + q);
p = cell(0);

%loop through LHC samples
for jj = 1:numsamples
    
    k = 1;                                      %init node counter
    %pick parameters within each range
    p{jj} = struct;                             %output parameter struct
    p{jj}.k1 = rng2lhs(k1rng,X(jj,1));          %TX rate
    p{jj}.k2 = rng2lhs(k2rng,X(jj,2));          %TL rate
    p{jj}.k3 = rng2lhs(k3rng,X(jj,3));          %degradation rate by protease
    p{jj}.delta1 = rng2lhs(delta1rng,X(jj,4));  %mRNA dilution/degradation
    p{jj}.delta2 = rng2lhs(delta2rng,X(jj,5));  %protein dilution
    p{jj}.RNAP = rng2lhs(RNAPrng,X(jj,6));      %total RNAP
    p{jj}.Ribo = rng2lhs(Riborng,X(jj,7));      %total ribosomes
    p{jj}.Ptot = rng2lhs(Ptotrng,X(jj,8));      %total protease
    p{jj}.DNA = rng2lhs(DNArng,X(jj,9:8+n));    %total DNA
    p{jj}.Kp = rng2lhs(Kprng,X(jj,9+n:8+2*n),'log');    %RNAP/DNA leaky binding cnst
    p{jj}.K2 = rng2lhs(K2rng,X(jj,9+2*n:8+3*n),'log');    %mRNA/ribosome binding cnst
    p{jj}.K3 = rng2lhs(K3rng,X(jj,9+3*n:8+4*n),'log');    %protein/protease binding cnst
    p{jj}.K0 = rng2lhs(K0rng,X(jj,9+4*n:8+4*n+m),'log');   %TF/DNA binding cnst
    p{jj}.K1 = rng2lhs(K1rng,X(jj,9+4*n+m:8+4*n+2*m),'log');   %RNAP/(TF/DNA) binding cnst
    p{jj}.s = rng2lhs(nrng,X(jj,9+4*n+2*m:8+4*n+3*m));          %coopertivity
    p{jj}.K02 = rng2lhs(K0rng,X(jj,9+4*n+2*m:8+4*n+3*m+q),'log');%TF/DNA binding cnst for two TFs at once
    
    %production parameters
    hparams = struct;
    %resource demand
    hparams.J = p{jj}.DNA./p{jj}.Kp.*(1+p{jj}.k1*p{jj}.RNAP./(p{jj}.K2*p{jj}.delta1));
    %leakiness
    hparams.T = p{jj}.k2*p{jj}.k1*p{jj}.RNAP*p{jj}.Ribo*p{jj}.DNA./(p{jj}.Kp.*p{jj}.K2.*p{jj}.delta1);
    %coopertivitity
    hparams.n = p{jj}.s;
    %loop through each edge
    for ii = 1:m
        if edge_weights(ii) > 0
            %activation; max: T*a/b; min: T
            hparams.a(ii) = max([p{jj}.Kp(k)/p{jj}.K1(ii),p{jj}.K1(ii)/p{jj}.Kp(k)])*...
                abs(edge_weights(ii))/p{jj}.K0(ii); %activation strength
            hparams.b(ii) = 1./p{jj}.K0(ii);        %1/binding constant
        else
            %repression; max: T; min: T*a/b
            hparams.a(ii) = min([p{jj}.Kp(k)/p{jj}.K1(ii),p{jj}.K1(ii)/p{jj}.Kp(k)])/...
                (abs(edge_weights(ii))*p{jj}.K0(ii)); %activation strength
            hparams.b(ii) = 1./p{jj}.K0(ii);        %1/binding constant
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
                    d = 1./(p{jj}.K02(s)*p{jj}.K0(rinds(r,1))) + 1./(p{jj}.K02(s)*p{jj}.K0(rinds(r,2)));
                    hparams.c(s) = max([p{jj}.Kp(k)/p{jj}.K1(ii),p{jj}.K1(ii)/p{jj}.Kp(k)])*d;
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
    gparams.k = p{jj}.k3;       %catalytic rate of degradation
    gparams.Ptot = p{jj}.Ptot;  %total amount of protease
    gparams.K = p{jj}.K3;       %binding constant for each node
    
    %combining into output parameter struct
    p{jj}.delta = p{jj}.delta2;     %dilution rate constant
    p{jj}.n = n;                %number of nodes
    p{jj}.m = m;                %number of edges
    p{jj}.h = hparams;          %production parameter struct
    p{jj}.a = hparams;          %production resource sharing struct
    p{jj}.g = gparams;          %degradation resource sharing struct
    
end

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


function out = rng2lhs(p_rng,lhsweight,scale)
%use weighting from latin hypercube to select parameters

if nargin < 3
    scale = '';
end

%create sample
if strcmp(scale,'log')
    out = 10.^(p_rng(1) + (p_rng(2) - p_rng(1)).*lhsweight(:)');
else
    out = p_rng(1) + (p_rng(2) - p_rng(1)).*lhsweight(:)';
end
