function [dz,Jout] = dynamicsbio_mulambda(t,z,funs,u,tfinal)

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
