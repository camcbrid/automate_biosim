function [p2,prng] = params_RS_LHS(A,numsamples)
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
X = lhsdesign(numsamples, 16+5*n+3*m+q);
p = cell(0);
p2 = cell(0);
%DNA, Kp K2 K3 x2
%loop through LHC samples
%can transform with icdf('dist_name',x)
for jj = 1:numsamples
    
    %NEED TO select parameters for resource sensor
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
    p{jj}.DNA = rng2lhs(DNArng,X(jj,9:10+n));    %total DNA
    p{jj}.Kp = rng2lhs(Kprng,X(jj,11+n:12+2*n),'log');    %RNAP/DNA leaky binding cnst
    p{jj}.K2 = rng2lhs(K2rng,X(jj,13+2*n:14+3*n),'log');    %mRNA/ribosome binding cnst
    p{jj}.K3 = rng2lhs(K3rng,X(jj,15+3*n:16+4*n),'log');    %protein/protease binding cnst
    p{jj}.K0 = rng2lhs(K0rng,X(jj,16+4*n:15+4*n+m),'log');   %TF/DNA binding cnst
    p{jj}.K1 = rng2lhs(K1rng,X(jj,16+4*n+m:15+4*n+2*m),'log');   %RNAP/(TF/DNA) binding cnst
    p{jj}.s = rng2lhs(nrng,X(jj,16+4*n+2*m:15+4*n+3*m));          %coopertivity
    p{jj}.K02 = rng2lhs(K0rng,X(jj,16+4*n+2*m:15+4*n+3*m+q),'log');%TF/DNA binding cnst for two TFs at once
    
    %parameters for each experiment
    p2{jj,1} = exp_split_paras(p{jj},1);
    p2{jj,2} = exp_split_paras(p{jj},2);
    p2{jj,3} = exp_split_paras(p{jj},3);
    p2{jj,4} = exp_split_paras(p{jj},4);
    p2{jj,5} = exp_split_paras(p{jj},5);
    
    %loop through each experiment
    for l = 1:5
        k = 1;                             	%init node counter
        %production parameters
        hparams = struct;
        %resource demand
        hparams.J = p2{jj,l}.DNA./p2{jj,l}.Kp.*(1+p2{jj,l}.k1*p2{jj,l}.RNAP./...
            (p2{jj,l}.K2*p2{jj,l}.delta1));
        %leakiness
        hparams.T = p2{jj,l}.k2*p2{jj,l}.k1*p2{jj,l}.RNAP*p2{jj,l}.Ribo*p2{jj,l}.DNA./...
            (p2{jj,l}.Kp.*p2{jj,l}.K2.*p2{jj,l}.delta1);
        %coopertivitity
        hparams.s = p2{jj,l}.s;
        
        %loop through each edge
        for ii = 1:m
            if all(isfield(p2{jj,l},{'Kp','K1','K0'}))
                if edge_weights(ii) > 0
                    %activation; max: T*a/b; min: T
                    hparams.a(ii) = max([p2{jj,l}.Kp(k)/p2{jj,l}.K1(ii),...
                        p2{jj,l}.K1(ii)/p2{jj,l}.Kp(k)])*...
                        abs(edge_weights(ii))/p2{jj,l}.K0(ii); %activation strength
                    hparams.b(ii) = 1./p2{jj,l}.K0(ii);        %1/binding constant
                else
                    %repression; max: T; min: T*a/b
                    hparams.a(ii) = min([p2{jj,l}.Kp(k)/p2{jj,l}.K1(ii),...
                        p2{jj,l}.K1(ii)/p2{jj,l}.Kp(k)])/...
                        (abs(edge_weights(ii))*p2{jj,l}.K0(ii)); %activation strength
                    hparams.b(ii) = 1./p2{jj,l}.K0(ii);        %1/binding constant
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
                            d = 1./(p2{jj,l}.K02(s)*p2{jj,l}.K0(rinds(r,1))) + ...
                                1./(p2{jj,l}.K02(s)*p2{jj,l}.K0(rinds(r,2)));
                            hparams.c(s) = max([p2{jj,l}.Kp(k)/p2{jj,l}.K1(ii),...
                                p2{jj,l}.K1(ii)/p2{jj,l}.Kp(k)])*d;
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
        gparams.k = p{jj}.k3;       %catalytic rate of degradation
        gparams.Ptot = p{jj}.Ptot;  %total amount of protease
        gparams.K = p{jj}.K3;       %binding constant for each node
        
        %combining into output parameter struct
        p2{jj,l}.delta = p{jj}.delta2;     %dilution rate constant
        p2{jj,l}.n = n;                %number of nodes
        p2{jj,l}.m = m;                %number of edges
        p2{jj,l}.h = hparams;          %production parameter struct
        p2{jj,l}.a = hparams;          %production resource sharing struct
        p2{jj,l}.g = gparams;          %degradation resource sharing struct
    end
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
    out = 10.^(min(p_rng) + (max(p_rng) - min(p_rng)).*lhsweight(:)');
else
    out = min(p_rng) + (max(p_rng) - min(p_rng)).*lhsweight(:)';
end


function p2 = exp_split_paras(p,expnum)

%NOTE Randomness in above parameters across experiments is bad. FIX

p2.k1 = p.k1;           %TX rate
p2.k2 = p.k2;           %TL rate
p2.k3 = p.k3;           %degradation rate by protease
p2.delta1 = p.delta1;   %mRNA dilution/degradation
p2.delta2 = p.delta2;   %protein dilution
p2.RNAP = p.RNAP;       %total RNAP
p2.Ribo = p.Ribo;       %total ribosomes
p2.Ptot = p.Ptot;       %total protease
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
        %module with constitutie RFP
        p2.DNA = p.DNA(1:end-1);    %total DNA
        p2.K0 = p.K0;               %TF/DNA binding cnst
        p2.K02 = p.K02;             %TF/DNA binding cnst for two TFs at once
        p2.K1 = p.K1;               %RNAP/(TF/DNA) binding cnst
        p2.Kp = p.Kp(1:end-1);      %RNAP/DNA leaky binding cnst
        p2.K2 = p.K2(1:end-1);      %mRNA/ribosome binding cnst
        p2.s = p.s;                 %coopertivity
        p2.K3 = p.K3(1:end-1);      %protein/protease binding cnst
end
