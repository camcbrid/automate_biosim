function A = combine_adjacency(A1,A2,varargin)
%combine adjacency matricies into one for simulating multiple modules

if nargin > 2
    A2 = combine_adjacency(A2,varargin{:});
end
if size(A1,1) < size(A1,2) || size(A2,1) < size(A2,2)
    error('Adjacency matrix wrong shape')
end

n1 = min(size(A1));
n2 = min(size(A2));
m1 = max(size(A1));
m2 = max(size(A2));

if n1 == m1 && n2 == m2
    %both square
    A = [A1(1:n1,1:n1),zeros(n1,n2);
        zeros(n2,n1),A2(1:n2,1:n2)];
elseif n1 ~= m1 && n2 ~= m2
    %both nonsquare
    A = [A1(1:n1,1:n1),zeros(n1,n2);
        zeros(n2,n1),A2(1:n2,1:n2);
        A1(n1+1:end,:),A2(n2+1:end,:)];
elseif n1 ~= m1 && n2 == m2
    %A1 nonsquare, A2 square
    A = [A1(1:n1,1:n1),zeros(n1,n2);
        zeros(n2,n1),A2(1:n2,1:n2);
        A1(n1+1:end,:),zeros(m1-n1,n2)];
else
    %A1 square, A2 nonsquare
    A = [A1(1:n1,1:n1),zeros(n1,n2);
        zeros(n2,n1),A2(1:n2,1:n2);
        zeros(m2-n2,n1),A2(n2+1:end,:)];
end
    