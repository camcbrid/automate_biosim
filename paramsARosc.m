function p = paramsARosc(A)
%assign parameters for production to each edge and degradation parameters
%for each protein. Takes in the weighted adacjency matrix, outputs
%parameter struct for functions h, g, a, and L to be used by makefuns.m

%defaults
if ~exist('A','var')
    error('not enough inputs')
end

%setup
B = (A ~= 0);               %logical adjacency matrix
n = min(size(A));           %number of nodes
m = sum(B(A ~= 0));         %number of edges

%output parameter struct
p = struct;
p.k3 = 270;             %degradation rate by protease
p.delta1 = 1;           %mRNA dilution/degradation
p.delta2 = [4,.5];      %protein dilution
p.Ptot = 0;             %total protease
p.K3 = [1200,1200];     %protein/protease binding cnst (n)
p.s = [2,4,2,4];        %coopertivity (m)

%production parameters
hparams = struct;
%resource demand
hparams.J = [0,0];
%leakiness
hparams.T = [20,2];
%coopertivitity
hparams.n = p.s;
%loop through each edge
hparams.a = [10,1,10,1];
hparams.b = [1,10,1,10];

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
