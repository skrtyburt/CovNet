function [pval, z_obs, z1, z2] = fisher_r2z_compare(r1, n1, r2, n2)

narginchk(4,4)

% Testing statistical significance of two correlations:

% Zobserved = (z1 – z2) / (square root of [ (1 / N1 – 3) + (1 / N2 – 3) ]

%convert r to z
z1 = 0.5*log((1+r1)/(1-r1));
z2 = 0.5*log((1+r2)/(1-r2));

% sample size norm
szab = sqrt(1/(n1-3) + 1/(n2-3));

% Z observed
z_obs = abs(z1-z2)/szab;

% significance testng (zcrit)
p = 1-normcdf(z_obs, 0, 1);
% two tailed p-value
	pval = 2*p;