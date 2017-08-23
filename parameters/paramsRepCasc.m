function p = paramsRepCasc(A)
%assign parameters for production to each edge and degradation parameters
%for each protein. Takes in the weighted adacjency matrix, outputs
%parameter struct for functions h, g, a, and L to be used by makefuns.m

%defaults
if ~exist('A','var')
    error('not enough inputs')
end

%output parameter struct
p = struct;
p.k1 = 100;             %TX rate
p.k2 = 350;             %TL rate
p.k3 = 270;             %degradation rate by protease
p.delta1 = 20;          %mRNA dilution/degradation
p.delta2 = 5;           %protein dilution
p.RNAP = 50;            %total RNAP
p.Ribo = 500;           %total ribosomes
p.Ptot = 0;             %total protease
p.DNA = [50,200];       %total DNA
p.K0 = [10,1];          %TF/DNA binding cnst
p.K1 = [1e9,1e9];       %RNAP/(TF/DNA) binding cnst
p.Kp = [2000,100];      %RNAP/DNA leaky binding cnst
p.K2 = [10000,1000];    %mRNA/ribosome binding cnst
p.K3 = [1200,1200];     %protein/protease binding cnst
p.s = [1,4];            %coopertivity


%setup
B = (A ~= 0);               %logical adjacency matrix
indvec = cumsum(sum(B,1));  %indicies of nodes in the edge counter
n = min(size(A));           %number of nodes
m = sum(B(A ~= 0));         %number of edges
k = 1;                      %init node counter
edge_weights = A(A ~= 0);   %extract weighting of each edge

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
    if any(ii == indvec)
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
