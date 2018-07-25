function [T,J] = calcTJ(p)
%takes in parameters for a circuit and returns vector of J's and T's with
%elements for each node in the circuit

%unpack parameters for circuit
RNAP = p.RNAP;          %total RNAP
Ribo = p.Ribo;          %total Ribosomes
k1 = p.k1;              %TX rate
k2 = p.k2;              %TL rate
delta1 = p.delta1;      %mRNA dilution
DNA = cell2mat(p.DNA);  %cell array of copy numbers
Kp = cell2mat(p.Kp);    %basal TX binding constant
K2 = cell2mat(p.K2);    %mRNA/ribo binding constant

%find J for each node in circuit
J = zeros(length(DNA),1);
T = zeros(length(DNA),1);
for ii = 1:length(DNA)
    J(ii) = DNA(ii)/Kp(ii)*(1+k1*RNAP/(delta1*K2(ii)));
    T(ii) = Ribo*RNAP*DNA(ii)*k1*k2/(Kp(ii)*K2(ii)*delta1);
end
