function J = runbio_jac(t,z,funs,u,mu,lambda)
%output Jacobian matrix for system with constant input, u

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
            ones(size(z))') + mu*kron(h(z,u),aJ(z,u)') + lambda*gJ(z) - L;
    else
        J = zeros(size(hJ(z,u)));
    end
else
    J = zeros(length(z));
end