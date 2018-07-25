function [F,Fprime] = Fhill(p,num_inputs)
%generate function handles for hill functions with different number of 
%inputs. Takes in parameters for one node and outputs Hill fcn handle and 
%function handle for the derivative

if nargin < 2
    num_inputs = 1;
    if nargin < 1
        p = struct('Kp',1000,'K0',1000,'K1',1000,'n',2);
    end
end

if num_inputs == 1
    %for single input hill function
    Kp = p.Kp;
    K0 = p.K0;
    K1 = p.K1;
    n = p.n;
    %[Kp,K0,K1,n] = deal(p{:});
    %Kp: basal dissociation constant RNAP with promoter
    %K1 : dissociation constant RNAP with TF/promoter complex (K)
    %K0 : dissociation constant of TF with promoter (k)
    a = Kp/(K1*K0);
    b = 1/K0;
    %output hill func and its derivative
    F = @(x) (1 + a*x.^n)./(1 + b*x.^n);
    Fprime = @(x) (n*(a-b)*x.^(n-1))./((1+b*x.^n).^2);
    
elseif num_inputs == 2
    %for two input hill function
    %input parameters
    Kp = p.Kp;
    K1 = p.K1;
    K0 = p.K0;
    n = p.n;
    [K11,K12] = deal(K1{:});
    [K01,K02,K012,K021] = deal(K0{:});
    [n1,n2] = deal(n{:});
    
    a1 = Kp/(K11*K01);
    a2 = Kp/(K12*K02);
    a3 = (Kp/K12)*(1/(K01*K012)+1/(K02*K021));
    b1 = 1/K01;
    b2 = 1/K02;
    b3 = (1/(K01*K012)+1/(K02*K021));
    
    num = @(x1,x2) (1 + a1*x1.^n1 + a2*x2.^n2 + a3*x1.^n1*x2.^n2);
    den = @(x1,x2) (1 + b1*x1.^n1 + b2*x2.^n2 + b3*x1.^n1*x2.^n2);
    
    %output hill function and both directional derivatives
    F = @(x1,x2) num(x1,x2)./den(x1,x2);
    Fprime1 = @(x1,x2) n1*x1.^(n1-1).*((a1+a3*x2.^n2).*(1+b2*x2.^n2)...
        -(b1+b3*x2.^n2).*(1+a2*x2.^n2))./den(x1,x2).^2;
    Fprime2 = @(x1,x2) n2*x2.^(n2-1).*((a2+a3*x1.^n1).*(1+b1*x1.^n1)...
        -(b2+b3*x1.^n1).*(1+a1*x1.^n1))./den(x1,x2).^2;
    Fprime = {Fprime1; Fprime2};
    
else
    %constant hill function
    F = @(x) 1;
    Fprime = @(x) 0;
end