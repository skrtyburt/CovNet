function [cov_mat, param_p_mat, perm_p_mat] = randshiftnull_cov(covTS,R)
narginchk(2,2)
% Input 
%   covTS - Regional PET values for a group as a 2D matrix where rows are
%           regions and columns are subjects.
%   R -     Number of permulations for the null distributions.
%
% Outputs
%   cov_mat -       emipirical covariance matrix (same as the one generated
%                   earlier in the workflow).
%   param_p_mat -   p-value matrix from the corr function.
%   perm_p_mat -    permutation p-value matrix.
%
% Version Control:
% Original version - Evgeny Chumin, Indiana University, 2023
%%
N = size(covTS,1); % row are regions/nodes
S = size(covTS,2); % columns are subjects

% generate the empirical matrix
[cov_mat, param_p_mat] = covariance(covTS,0);
% mask of positive correlations
pos_mat=cov_mat;
pos_mat(pos_mat<0)=0;
pos_mask = logical(pos_mat);
% mask of negative correlations
neg_mat=cov_mat;
neg_mat(neg_mat>0)=0;
neg_mask=logical(neg_mat);

% randshift null
% permute intensity values independently within each subject
pcnt = zeros(N,N);
for r=1:R
    tsr = zeros(N,S);
    for s=1:S
        tsr(:,s) = covTS(randperm(N),s);
    end
    % genrate the null matrix
    covCorr_r = covariance(tsr,0);
    % separate empirical data postive and negative edges 
    rpos = covCorr_r;
    rpos(~pos_mask)=0;
    rneg = covCorr_r;
    rneg(~neg_mask)=0;
    
    % idependently assess if emp+ > null+
    pcnt = pcnt+(pos_mat>rpos);
    % emp- < null-
    pcnt = pcnt+(neg_mat<rneg);
end

% pval
perm_p_mat = 1.-pcnt./R;
end