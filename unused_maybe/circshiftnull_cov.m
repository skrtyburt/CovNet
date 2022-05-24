% Adapted from a rotational shift null model implementation of Dr. Olaf
% Sporns, Indiana University, to operate on subject series in order to
% estimate a permutation p-value matrix for potnetial thresholding of the
% covariance matrix. 
function [cov_mat, param_p_mat, perm_p_mat, thr_cov_mat] = circshiftnull_cov(covTS,R,pthr)

N = size(covTS,2);
lts = size(covTS,1);

offsets = [-lts:1:lts];

[cov_mat, param_p_mat] = corr(covTS); % Pearson

% circshift null
pcnt = zeros(N,N);
for r=1:R
    tsr = zeros(lts,N);
    for n=1:N
        tsr(:,n) = circshift(covTS(:,n),offsets(randi(length(offsets))));
    end
    covCorr_r = corr(tsr,'type','spearman');
    
    pcnt = pcnt+(cov_mat>covCorr_r);       
end

% pval
perm_p_mat = 1.-pcnt./R;

p = perm_p_mat;
p(p>pthr)=0; p=logical(p);

thr_cov_mat = cov_mat;
thr_cov_mat(~p)=nan;