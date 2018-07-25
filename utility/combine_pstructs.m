function p = combine_pstructs(p1,p2)
%Combine parameter structs into one for combining modules in simulation.
%Arguments should be structs (can be nested structs). Takes cellular
%properties from first argument.

p = struct();
fields1 = fieldnames(p1);
fields2 = fieldnames(p2);
fields = unique([fields1(:);fields2(:)]);

%loop through all fields
for i = 1:length(fields)
    if ~isfield(p1,fields{i})
        p1.(fields{i}) = [];
    end
    if ~isfield(p2,fields{i})
        p2.(fields{i}) = [];
    end
    %if parameter is a cell property, take values from first parameter
    %vector
    if strcmp(fields{i},'k1') || strcmp(fields{i},'k2') || strcmp(fields{i},'k3')...
            || strcmp(fields{i},'delta1') || strcmp(fields{i},'delta2') || ...
            strcmp(fields{i},'delta') || strcmp(fields{i},'RNAP') || ...
            strcmp(fields{i},'Ribo') || strcmp(fields{i},'Ptot') || strcmp(fields{i},'k')
        p.(fields{i}) = p1.(fields{i});
    else
        if ~isnumeric(p1.(fields{i})) || ~isnumeric(p2.(fields{i}))
            %check if p.(field{i}) is a struct, and, if so, use recursion
            p.(fields{i}) = combine_pstructs(p1.(fields{i}),p2.(fields{i}));
        else
            if strcmp(fields{i},'n') || strcmp(fields{i},'m') || strcmp(fields{i},'m2')
                %combine numeric arrays
                p.(fields{i}) = p1.(fields{i})(:) + p2.(fields{i})(:);
            else
                %combine numeric arrays
                p.(fields{i}) = [p1.(fields{i})(:); p2.(fields{i})(:)]';
            end
        end
    end
end
