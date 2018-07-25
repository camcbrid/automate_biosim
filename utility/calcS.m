function s = calcS(y, u, A)
%A = calcS(Y, U, A)
%calculate the derivatives of each node with respect to the disturbance J_0
%h0 represents change in node from direct resource perturbation 
%g represents change to node from change in input concentration
%hij represents change to node from change in total concentration 
%/avaialable resource
%
%Y is the 


%adjency matrix (for a cascade)
%A = diag(ones(2,1),1);
%u = logspace(0,10,25)';

n = min(size(A));

%input vectors
w = cell(0);
w{1} = y.x{1}.y5(:,[1,1:2])*A;
w{1}(:,1) = u;
w{2} = y.x{2}.y5(:,[1,1:2])*A;
w{2}(:,1) = u;
w{3} = y.x{3}.y5(:,[1,1:2])*A;
w{3}(:,1) = u;

%parameters
T = y.T;
J = y.J;
delta2 = y.p.delta2;
Fp = y.Fp;
F = y.F;
Q = y.Qcalc;
J0 = y.J0;

xss = cell(0);      denom = cell(0);
g = cell(0);        h0 = cell(0);       h = cell(0);
G = cell(0);        H = cell(0);
s = cell(0);
Delta = cell(0);    f = cell(0);

%loop through modules
for k = 1:3
    %Q{k} = J(1)*F1(uvec) + J(2)*F2(xA.y4(:,1)) + J(3)*F3(xA.y4(:,2));
    denom{k} = (1+Q{k}+J0{k});
    %loop through nodes in module
    for ii = 1:n
        u2 = w{k}(:,ii);
        xss{k}(:,ii) = T{k}(ii)*F{k,ii}(u2)./(delta2*denom{k});
        g{k}(:,ii) = xss{k}(:,ii).*Fp{k,ii}(u2)./F{k,ii}(u2);
        h0{k}(:,ii) = xss{k}(:,ii)./denom{k};
        %loop through nodes in module
        for jj = 1:n
            u3 = w{k}(:,jj);
            h{k}(:,ii,jj) = (xss{k}(:,ii).*J{k}(jj).*Fp{k,jj}(u3))./denom{k};
        end
    end
    %HELP
    H{k} = shiftdim(h{k},1);
    %loop through levels of induction
    for m = 1:length(u)
        G{k}(:,:,m) = diag(g{k}(m,:));
        %calculate du_j/dJ0
        s{k}(:,m) = -(eye(n) + (H{k}(:,:,m) - G{k}(:,:,m))*A')\...
            h0{k}(m,:)';
        %HELP
        %for l = 1:n
        %    f{k,l} = J{k}(l)*Fp{k,l}();
        %    Delta{k} = f{k,l}'*s{k};
        %end
    end
end
