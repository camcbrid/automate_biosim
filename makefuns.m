function [funs,p] = makefuns(A,parafun)
%make function handles for production, resource sharing, and degradation
%dynamics. Can be used in runbio.m run by ode15s to run dynamics of a
%protein system with QSS assumption. Inputs are adjacency matrix and
%parameter function that output parameter struct and takes adjacency matrix
%as input. Also returns parameter struct used for debugging

%defaults
if nargin < 2
    parafun = @params_dist;
    if nargin < 1
        error('not enough inputs')
    end
end

%create parameter struct, either randomly or set individual parameters
p = parafun(A);
%number of nodes
n = min(size(A));
%dilution matrix
if length(p.delta) == 1
    L = diag(p.delta*ones(n,1));
else
    L = diag(p.delta);
end
%function handles for dynamics
funs = struct;
funs.h = @(x,u) hfun(x,p.h,u,A);        %production function
funs.a = @(x,u) afun(x,p.a,u,A);        %production resource sharing
funs.g = @(x) gfun(x,p.g);              %degradation including resource sharing
funs.L = L;                             %dilution matrix
funs.f = @(x,u,mu,lambda) hfun(x,p.h,u,A).*((afun(x,p.a,u,A) - ones(size(x)))*mu + ...
    ones(size(x))) + lambda*gfun(x,p.g) - L*x;  %complete dynamics
%Jacobians
funs.hJ = @(x,u) hJfun(x,p.h,u,A);      %production function Jacobian
funs.aJ = @(x,u) aJfun(x,p.a,u,A);      %production resource sharing Jacobian (vector)
funs.gJ = @(x) gJfun(x,p.g,A);          %degradation Jacobian (complete clique)
funs.J = @(x,u,mu,lambda) hJfun(x,p.h,u,A).*kron((afun(x,p.a,u,A) - ...
    ones(length(x),1))*mu + ones(length(x),1),ones(size(x))') + ...
    mu*kron(hfun(x,p.h,u,A),aJfun(x,p.a,u,A)') + lambda*gfun(x,p.g) - L; %Jacobian
%helper functions
funs.F = @(x,a,b,n) F2(x,a,b,n);                %Hill function
funs.Fp = @(x,a,b,s,n,inds) Fp2(x,a,b,s,n,inds);%Hill function derivative w/r/t all states
funs.makeU = @(x,u,A) makeU(x,A,u);             %create inputs for each node
funs.combineFU = @(U,p,A) combineFU(U,p,A);     %evalulate Hill functions at U's
funs.combineFpU = @(U,p,A,inds) combineFpU(U,p,A,inds); %eval Hill fcn derivatives



function out = hfun(x,p,u,A)
%make h function representing production and accounting for resource
%sharing of RNAP and ribosomes

%make adjacency matrix all ones and zeros
A = (A ~= 0);
%make input array, U
U = makeU(x,A,u);
%evalulate normalized Hill functions for each input
Fout = combineFU(U,p,A);
%combine T's and F's
out = cell2mat(cellfun(@(a,b) a.*b, num2cell(p.T(:)), Fout(:),...
    'UniformOutput',false));



function out = afun(x,p,u,A)
%make alpha function for resource sharing terms

%make adjacency matrix all ones and zeros
A = (A ~= 0);
%make input array, U
U = makeU(x,A,u);
%evalulate normalized Hill functions for each input
Fout = combineFU(U,p,A);
%combine F's and J's
normvec = cell2mat(cellfun(@(x,y) x.*y, num2cell(p.J(:)), Fout(:),...
    'UniformOutput',false));
%alpha (can be scalar or vector)
out = 1./(1+sum(normvec));



function out = gfun(x,p)
%make degradation function. Outputs a vector with same size as x. x is the
%state vector, p is the parameter struct for degradation, and A is the
%adacency matrix. Degradation applies to all proteins with this model
%x input must be row vec of col vecs [x0,x1,x2,...]

%parameters
k = p.k;
P = p.Ptot;
K = p.K;
%make g function
xnorm = zeros(size(x));
for ii = 1:size(x,1)
    xnorm(ii,:) = x(ii,:)/K(ii);    %protein concentration over their MM cnsts
end
out = -k*P*xnorm./(1 + sum(xnorm,1));



function [U,inds] = makeU(x,A,u)
%creates the input array, U. x is the state vector, A is the adacency
%matrix, and u is the external input
if nargin < 3
    u = 0;
    error('need external input, u')
end
%init
inds = cell(0);
n = min(size(A));       %number of nodes
U = cell(n,1);
%loop through inputs for each node to create input cell variable, U
%U{i} contains the concentration of the input to node {i}
mask = 1:size(A,1);
for q = 1:n
    indvec = mask.*A(:,q)';           %indicies of inputs to node q
    inds{q} = indvec(indvec ~= 0);    %eliminate 0's
    if any(inds{q} > n)
        %if node has external input (ind > n), then U is input states and
        %external input. External input, u, appears in Hill fcns like other
        %inputs
        if ~isempty(inds{q}(inds{q}<=n))
            U{q} = [x(inds{q}(inds{q}<=n),:);u*ones(1,size(x,2))];
        else
            U{q} = u*ones(1,size(x,2));
        end
    elseif ~isempty(x(inds{q}))
        %input is state vector based off connection
        U{q} = x(inds{q},:);
    else
        %no input to node q (constitutive)
        U{q} = num2cell(zeros(1,size(x,2)));
    end
end



function Fout = combineFU(U,p,A)
%creates Hill function handles and evalulates each at the input point, U. U
%is the input set (cell array of states), p is the parameter struct for
%production, and A is the adjacency matrix.

%parameters
a = p.a;            %numerator parameters
b = p.b;            %denominator parameters
s = p.n;            %coopertivity
n = min(size(A));   %number of nodes
invec = sum(A,1)';  %out degree
%init
k = 1;              %counter for edge parameters
F = cell(0);
Fout = cell(0);
for ii = 1:n
    if invec(ii) >= 1
        %n inputs to node i
        F{ii} = @(x) F2(x,a(k:k+invec(ii)-1),b(k:k+invec(ii)-1),s(k:k+invec(ii)-1));
        k = k + invec(ii);
    else
        %constant output if no input
        F{ii} = @(x) ones(size(x,2));
    end
    %evalulate each F at the input
    Fout{ii} = F{ii}(U{ii});
end



function out = F2(x,a,b,n)
%make generic Hill function ignoring combinatoral promoter binding.
%x is the input vector, a is the Hill function numerator parameter vector
%(cell array), b is the Hill function denominator parameter vector (cell
%array), s is the coopertivity vector (cell). Outputs a scalar

%init
numvec = zeros(size(x));
denvec = zeros(size(x));
%loop through each input,x
for ii = 1:size(x,1)
    numvec(ii,:) = a(ii)*x(ii,:).^n(ii);
    denvec(ii,:) = b(ii)*x(ii,:).^n(ii);
end
%Hill function output
out = (1+sum(numvec,1))./(1+sum(denvec,1));



function out = hJfun(x,p,u,A)
%make Jacobian dh/dx for nominal protein production without resources
%sharing

%make adjacency matrix all ones and zeros
A = (A ~= 0);
%make input array, U
[U,inds] = makeU(x,A,u);
%derivatives of Hill fcns, F (outputs a matrix)
Fpout = combineFpU(U,p,A,inds);
%multiply each row by the node's T
out = p.T(:).*Fpout;



function out = aJfun(x,p,u,A)
%make Jacobian da/dx for resource sharing terms (vector)

%make adjacency matrix all ones and zeros
A = (A ~= 0);
%make input array, U
[U,inds] = makeU(x,A,u);
%Hill functions and derivatives with inputs
Fout = combineFU(U,p,A);
Fpout = combineFpU(U,p,A,inds);
%intermediate
normvec = cell2mat(cellfun(@(x,y) x.*y, num2cell(p.J(:)), Fout(:),...
    'UniformOutput',false));
%Jacobian of a(x)
out = - sum(p.J(:).*Fpout,2)/(1 + sum(normvec,1)).^2;



function out = gJfun(x,p,A)
%make the Jacobian dg/dx of the degradation function

A = (A ~= 0);       %make adjacency matrix all ones and zeros
n = min(size(A));   %number of nodes
%parameters
k = p.k;
Ptot = p.Ptot;
K = p.K;
%init
xnorm = zeros(n,1);
Knorm = zeros(n,1);
%loop through each node, creating intermediates
for ii = 1:n
    xnorm(ii) = x(ii)/K(ii);
    Knorm(ii) = 1/K(ii);
end
den = 1 + sum(xnorm);               %denominator
outvec = den*Knorm;                 %diagonal entries...
B = diag(outvec);                   %into a diagonal matrix
C = kron(xnorm,Knorm');             %rest of entries (complete)
out = -k*Ptot*(B - C)./(den.^2);    %dg/dx



function out = Fp2(x,a,b,s,n,inds)
%make Hill function derivative output vector. Returns a vector of the
%derivative with respect to each input. x is the input state vector, a is
%the Hill function numerator parameter vector (cell array), b is the Hill
%function denominator parameter vector (cell array), s is the coopertivity
%vector (cell), n is the number of nodes in the network, and inds is the
%indicies of the inputs. If indicies is shorter than x, then some x's are
%external inputs.

%init
numvec = zeros(length(x),1);
denvec = zeros(length(x),1);
vals = zeros(length(inds),1);
out = zeros(1,n);               %output a row of length n
%loop through each input,x
for ii = 1:min([length(a),length(x),length(b),length(s)])
    numvec(ii) = a(ii)*x(ii).^s(ii);
    denvec(ii) = b(ii)*x(ii).^s(ii);
end
asum = (1+sum(numvec));
bsum = (1+sum(denvec));
%calculate the derivative w/r/t each index at one point, x
for jj = 1:length(inds)
    vals(jj) = s(jj)*x(jj)^(s(jj)-1)*(a(jj)*bsum - b(jj)*asum)/(bsum^2);
end
%put vals into correct spots in row in Jacobian
out(inds) = vals;



function Fpout = combineFpU(U,p,A,inds)
%creates derivative of Hill function handles and evalulates each at the
%input point, U. U is the input set (cell array of states), p is the
%parameter struct for production, A is the adjacency matrix, and inds is
%the cell array holding the indicies that input to node ii.

%init
n = min(size(A));           %number of nodes
invec = sum(A,1)';          %out degree
k = 1;                      %counter for edges
Fpout = zeros(n);           %Jacobian dh/dx
inds2 = cell(0);            %index vector without external inputs
%parameters
a = p.a;
b = p.b;
s = p.n;
%derivative of hill functions
Fp = cell(0);
%Hill functions
for ii = 1:n
    inds2{ii} = inds{ii}(inds{ii} <= n); %chop off external inputs indicies
    if ~isempty(inds2{ii})
        %if there are some inputs from other states
        numext = sum(inds{ii} > n); %number of external inputs to node ii
        %row in the Jacobian
        Fp{ii} = @(x) Fp2(x,a(k:k+invec(ii)-1-numext),b(k:k+invec(ii)-1-...
            numext),s(k:k+invec(ii)-1-numext),n,inds2{ii});
    else
        %only input is external to node ii
        Fp{ii} = @(x) zeros(1,n);
    end
    k = k + invec(ii);  %advance edge counter
    %build Jacobian dh/dx one row at a time
    Fpout(ii,:) = Fp{ii}(U{ii});
end

