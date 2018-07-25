function [y, grad] = dynamicsbio_fmincon(t, z, funs, u, mu, lambda)
%[Y, GRAD] = dynamicsbio_fmincon(T, Z, FUNS, U, MU, LAMBDA)
%Simulate bio process with Hill function kinetics (simplified) with
%constant input, u for solving for the equilibrium point(s) using fmincon.

%T is the time input, Z is the state vector, FUNS is the
%struct of function handles returned by makefuns.m, U is the external
%input, MU is a parameter \in [0,1] controlling production resource
%sharing, and LAMBDA is a parameter \geq 0 controlling degradation resource
%sharing. Y is the sum squared value of the derivative vector DZ and
%GRAD is the gradient of Y with respect to the state vector Z.

%defaults
if nargin < 5
    mu = 1;             %turn on production resource sharing [0,1]
    if nargin < 6
        lambda = 1;     %turn on degradation resource sharing [0,1]
    end
end

%dynamics functions
h = funs.h;     %production
a = funs.a;     %production resource sharing
g = funs.g;     %degradation with resource sharing
L = funs.L;     %dilution matrix

%dynamics
dz = h(z,u).*((a(z,u) - ones(size(z)))*mu + ones(size(z))) + lambda*g(z) - L*z;
%objective function scalar output
y = rssq(dz).^2;

%Jacobian
if nargout == 2 && all(isfield(funs,{'hJ','aJ','gJ'}))
    %Jacobians for each function, h, a, and g
    %vectorizing doesn't work here
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
    grad = 2*sum(dz(:).*J,1);
end
