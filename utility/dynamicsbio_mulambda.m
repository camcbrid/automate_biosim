function [dz, Jout] = dynamicsbio_mulambda(t, z, funs, u, tfinal)
%[DZ, JOUT] = dynamicsbio_mulambda(T, Z, FUNS, U, TFINAL)
%Returns derivative for bio system with dynamics given by FUNS by slowly
%sweeping from (MU,LAMBDA) = (0,0) to (MU,LAMBDA) = (1,1).
%
%T is the time input, Z is the state vector, FUNS is the struct of function
%handles returned by makefuns.m, U is the external input, MU is a parameter
%\in [0,1] controlling production resource sharing, and LAMBDA is a
%parameter \geq 0 controlling degradation resource sharing. DZ is the
%time derivative of the system evalulated at the state Z. J is the Jacobian
%of the system evalulated at the state Z. If the fields hJ, gJ, and aJ are
%not present in FUNS or if FUNS.a is non-scalar, J is returned as a matrix
%of zeros.

%dynamics functions
h = funs.h;     %production
a = funs.a;     %production resource sharing
g = funs.g;     %degradation with resource sharing
L = funs.L;     %dilution matrix

mu = t/tfinal;         %turn on production resource sharing [0,1]
lambda = t/tfinal;     %turn on degradation resource sharing [0,1]

%dynamics
dz = h(z,u).*((a(z,u)-ones(size(z)))*mu+ones(size(z))) + lambda*g(z) - L*z;

%Jacobian
if nargout == 2 && all(isfield(funs,{'hJ','aJ','gJ'}))
    %Jacobians for each function, h, a, and g
    hJ = funs.hJ;
    gJ = funs.gJ;
    aJ = funs.aJ;
    
    %works for vector version of da/dx (not matrix)
    if min(size(aJ(z,u))) == 1
        Jout = hJ(z,u).*kron((a(z,u)-ones(length(z),1))*mu + ones(length(z),1),...
            ones(size(z))') + mu*kron(h(z,u),aJ(z,u)) + lambda*gJ(z) - L;
    else
        Jout = zeros(size(hJ(z,u)));
    end
end
