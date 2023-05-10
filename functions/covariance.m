function [cov_mat, cov_pval_mat] = covariance(data,flag,roi_labels,lims,cmap,covs)

% data -        regions(rows) by subject(column)
% covs -        covariates to account for in subject-by-covariate matrix
% roi_labels -  a cell array of string lables for regions/nodes
% flag -        if set to 1, generates a figure of the covariance patrix
% lims -        [# #] setting the lower and upper bounds of the figure
% cmap -        colormap

zdata = zscore(data)';

if exist('covs','var')
    covs = covs';
    [cov_mat, cov_pval_mat] = partialcorr(zdata,covs,'type','pearson');
else
    [cov_mat, cov_pval_mat] = corr(zdata,'type','pearson');
end

N = length(cov_mat);

if exist('flag','var') && flag == 1
% show the fully connected correlation adjacency matrix
    figure(...
        'units','inches',...
        'position',[1 1 6 7.5],...
        'paperpositionmode','auto');
    subplot(3,1,1)
    imagesc(zdata)
    ylabel('Subjects'); xlabel('Regions')
    title('Z-Scored Series')
    colorbar
    subplot(3,1,2:3)
    imagesc(cov_mat); axis square
    if exist('cmap','var')
        colormap(cmap);
    end
    if exist('lims','var')
        caxis(lims)
    end
    if exist('roi_labels','var') && length(roi_labels)==N
        xticks(1:1:N); yticks(1:1:N);
        yticklabels(roi_labels); xticklabels(roi_labels); xtickangle(45)
    end
    colorbar
    title('Pearson Covariance Matrix')
end