function [cov_mat, param_p_mat, perm_p_mat] = randshiftnull_cov(covTS,R,covs)
narginchk(2,3)
N = size(covTS,1); % row are regions/nodes
S = size(covTS,2); % columns are subjects

if exist('covs','var')
    [cov_mat, param_p_mat] = covariance(covTS,0,covs);
else
    [cov_mat, param_p_mat] = covariance(covTS,0);
end
pos_mat=cov_mat;
pos_mat(pos_mat<0)=0;
pos_mask = logical(pos_mat);
neg_mat=cov_mat;
neg_mat(neg_mat>0)=0;
neg_mask=logical(neg_mat);

% randshift null
pcnt = zeros(N,N);
for r=1:R
    tsr = zeros(N,S);
    for s=1:S
        tsr(:,s) = covTS(randperm(N),s);
    end
    if exist('covs','var')
        covCorr_r = covariance(tsr,0,covs);
    else
        covCorr_r = covariance(tsr,0);
    end
    rpos = covCorr_r;
    rpos(~pos_mask)=0;
    rneg = covCorr_r;
    rneg(~neg_mask)=0;
            
    pcnt = pcnt+(pos_mat>rpos);
    pcnt = pcnt+(neg_mat<rneg);
end

% pval
perm_p_mat = 1.-pcnt./R;
end