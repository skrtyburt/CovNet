function [perm_pmat] = perm_ttest2_cov(data1, data2,R)

narginchk(3,3)

[N1,S1] = size(data1);
[N2,S2] = size(data2);

if N1 ~= N2
    fprintf(2,'Unequal number of regions on rows\n')
    return
end

gvec = ones(1,S1);
gvec(1,S1+1:S1+S2)=2;
S = length(gvec);

dataall = horzcat(data1,data2);

C1 = covariance(data1);
C2 = covariance(data2);
Ediff = abs(C1-C2);

pcnt = zeros(N1,N1);

for r = 1:R
    gperm = gvec(randperm(S));
    P1 = covariance(dataall(:,gperm==1));
    P2 = covariance(dataall(:,gperm==2));
    Pdiff = abs(P1-P2);
    pcnt = pcnt + (Ediff>Pdiff);
    clear P1 P2 gperm
end

perm_pmat = 1 - (pcnt./R);
