function [perm_pmat, mxcor_pmat] = perm_ttest2_cov(data1, data2,R)
%%    Permulations testing of edge-level differences    %%
% Input 
%   data1 - 2D Region-by-Subject matrix for group1
%   data2 - 2D Region-by-Subject matrix for group2
%   R -     Number of group affiliation permutations
%
% Output
%   perm_pmat -     Matrix of permutation p-values.
%   mxcor_pmat -   Maximum statistic corrected p-values
%
% Version Control:
% Original version - Evgeny Chumin, Indiana University, 2022
%%
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
mxdiff = zeros(N1,N1);

for r = 1:R
    gperm = gvec(randperm(S));
    P1 = covariance(dataall(:,gperm==1));
    P2 = covariance(dataall(:,gperm==2));
    Pdiff = abs(P1-P2);
    mxdiff=max(cat(3,mxdiff,Pdiff),[],3);
    pcnt = pcnt + (Ediff>Pdiff);
    clear P1 P2 gperm
end
perm_pmat = 1 - (pcnt./R);

m= triu(logical(ones(N1)),1);
mxdist = mxdiff(m);
n=length(mxdist);
mxcor_pmat = arrayfun(@(x) (1-sum(x>mxdist)/n),Ediff);



