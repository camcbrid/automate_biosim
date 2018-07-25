function [A, n] = augment_input(B)
%[A, n] = augment_input(B)
%add input to first node in adjacency matrix
%
%B is the input adjacency matrix or cell array of adjacency matricies. n is
%the number of nodes in the system, and A is the output adjacency matrix
%with an external input appended to the (n+1)st row of B. If B is a cell
%array, A and n outputs as cell array with the same number of elements as
%B.

if iscell(B)
    %apply to all arrays within the cell array
    A = cell(0);        %output adjacency matrix
    n = cell(0);        %number of nodes in B
    for i = 1:length(B)
        [A{i},n{i}] = augment_input(B{i});
    end
    return
end

n = min(size(B));                   %number of nodes
A = [B; 1, zeros(1,size(B,2)-1)];   %append [1,0,0,0,...] to last row
