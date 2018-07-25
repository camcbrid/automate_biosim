function J = dynamicsbio_jac(t, z, funs, u, mu, lambda)
%J = dynamicsbio_jac(T, Z, FUNS, U, MU, LAMBDA)
%Output Jacobian matrix for system with constant input, u. 

%T is the time input, Z is the state vector, FUNS is the struct of function
%handles returned by makefuns.m, U is the external input, MU is a parameter
%\in [0,1] controlling production resource sharing, and LAMBDA is a 
%parameter \geq 0 controlling degradation resource sharing. J is the 
%Jacobian of the system evalulated at the state Z. If the fields hJ, gJ, 
%and aJ are not present in FUNS or if funs.a is non-scalar, J is returned 
%as a matrix of zeros.

%defaults
if nargin < 5
    mu = 1;
    if nargin < 6
        lambda = 1;
    end
end

%Jacobian
if all(isfield(funs,{'hJ','aJ','gJ'}))
    %Jacobians for each function, h, a, and g
    h = funs.h;
    a = funs.a;
    g = funs.g;
    L = funs.L;
    hJ = funs.hJ;
    gJ = funs.gJ;
    aJ = funs.aJ;
    
    %works for vector version of da/dx (not matrix)
    if min(size(aJ(z,u))) == 1
        J = hJ(z,u).*kron((a(z,u)-ones(length(z),1))*mu + ones(length(z),1),...
            ones(size(z))') + mu*kron(h(z,u),aJ(z,u)) + lambda*gJ(z) - L;
    else
        J = zeros(size(hJ(z,u)));
    end
else
    J = zeros(length(z));
end