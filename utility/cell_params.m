function p = cell_params(p1, varargin)
%P = cell_params(P1, varargin)
%set cellular properties 'k1', 'k2', 'k3', 'RNAP', 'Ribo', 'Ptot',
%'delta1', 'delta2', and 'delta' to be the same for all parameter structs. 
%Take cell properties from first parameter struct.
%
%P1 is a struct containing parameters

if nargin > 1
    p1 = [{p1}, varargin(:)'];
end

%fields that should be the same across all parameter structs
cellfields = {'k1','k2','k3','RNAP','Ribo','Ptot','delta1','delta2','delta'};

%set cellfields to the same in all structs
p = p1;
for i = 1:length(p1)
    for j = 1:length(cellfields)
        p{i}.(cellfields{j}) = p1{1}.(cellfields{j});
    end
end