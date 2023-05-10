function covariance_analysis_tier1(cellData,ROIlabels,cmap,pval,covariates)
narginchk(4,5)
% Input 
%   cellData - cell variable where rows and columns are grouping variables
%              (i.e. age or sex) for a group (i.e. control or mutant).
%              -First row and column of cell should have string labels for
%              groups.
%              -Value {1,1} should contain the group label (i.e. Control).
%              -Data should have regions along rows and subjects/animals
%              along colums.
%
%   ROIlabels - a cell row of string labels for regions (must be the same
%               size as number of rows in each of the data cells).
%   
%   cmap - colormap.
%
%   pval - a single or range of p-values for edge level thresholding of
%   networks, based on permutation testing of significance of correlations
%
%   covariates - a cell structure same size as cellData with the same
%   labels in first row/column, where data are organized in the same order
%   as cell data: each row is a covaraite with subjects in columns ordered
%   the same as in cellData
%
% Output
%
%
% Version Control:
% Original version - Evgeny Chumin, Indiana University, 2022
%%
%% Perform data checks
% cell structure check variable names in first row col
 grp = cellData{1,1};       % group label
 rNames = cellData(:,1);    % row subgroup labels
 cNames = cellData(1,:);    % column subgroup labels
 
 % check that the data are region by subject by making sure number of rows
 % matched number of regions in labels variable.
 N = length(ROIlabels);     % number of regions
 % there is a +1 for (1,1) which is why most loops start at 2.
 Nr = length(rNames);       % number of row groupings
 Nc = length(cNames);       % number of column groupings
 
 % count number of regions and animals in each data cell
 for r=2:Nr % every row
     for c=2:Nc % every column
         [R(r-1,c-1), S(r-1,c-1)] = size(cellData{r,c});     % # of regions and subjects
     end
 end
 
if length(unique(R))>1   % if number of regions varies 
    fprintf(2,'Number of regions varies. Check that cellData was made correctly.\n')
    return % exit
end
if R(1,1) ~= length(ROIlabels) % if number of regions does not match nubmer of labels
    fprintf(2,'Number of labels does not match number of regions. Exiting...\n')
    return % exits
end

%% Generate covariance matrices
% Set labels
covMat = rNames;        covMat(1,1:Nc) = cNames;        % covariance matrices
cov_paramP = rNames;    cov_paramP(1,1:Nc) = cNames;    % parametric p value matrices

outroot = ['ca_tier1_figures_' grp];

ver = length(dir(fullfile(pwd,[outroot '*'])));
if ver == 0
    outdir = fullfile(pwd,outroot);
else
    outdir = fullfile(pwd,[outroot '_run' num2str(ver+1)]);
end
clear ver
if ~exist(outdir,'dir')
    mkdir(outdir)
end

for r=2:Nr % every row
    for c=2:Nc % every column
        % generate covariance matrices
        if exist('covariates','var')
            [covMat{r,c}, cov_paramP{r,c}]=covariance(cellData{r,c},1,ROIlabels,[-1 1],cmap,covariates{r,c});
            sgtitle([grp ' ' cellData{r,1} ' ' cellData{1,c} ' partialCorr'])
        else
            [covMat{r,c}, cov_paramP{r,c}]=covariance(cellData{r,c},1,ROIlabels,[-1 1],cmap);
            sgtitle([grp ' ' cellData{r,1} ' ' cellData{1,c}])
        end
        
        filename = fullfile(outdir, ['covMat_' grp '_' cellData{r,1} '_' cellData{1,c} '.pdf']);
        print(gcf,'-dpdf',filename)
        clear filename
    end
end
close all
% Parametric p-values of correlation significance are saved for
% completeness, but the recommended methods for identifying signifance edge
% covariance is permutations testing.

%% Permutation testing of significant edge covariance 
% 10,000 permutations
cov_permP = rNames;    cov_permP(1,1:Nc) = cNames;    % permutation p value matrices
for r=2:Nr % every row
    for c=2:Nc % every column
        if exist('covariates','var')
            [~,~,cov_permP{r,c}] = randshiftnull_cov(cellData{r,c},10000,covariates{r,c});
        else
            [~,~,cov_permP{r,c}] = randshiftnull_cov(cellData{r,c},10000);
        end
    end
end

%% community detection 
% mrcc - Jeub 2018
agrMat = rNames;    agrMat(1,1:Nc) = cNames;    % multiscale agreement matrices
mrccPartition = rNames;    mrccPartition(1,1:Nc) = cNames;    % multiresolution consensus modular parition
allPartitions = rNames;    allPartitions(1,1:Nc) = cNames;    % all hierarchical partitions
for r=2:Nr % every row
    for c=2:Nc % every column
        % multiresolution concenus clustering 
        [agrMat{r,c}, mrccPartition{r,c},allPartitions{r,c}] = mrcc_wrapper(covMat{r,c},10000,1,ROIlabels);
        mrccPartition{r,c} = fcn_relabel_partitions(mrccPartition{r,c});
        sgtitle([grp ' ' cellData{r,1} ' ' cellData{1,c}])
        filename = fullfile(outdir, ['mrcc_out_' grp '_' cellData{r,1} '_' cellData{1,c} '.pdf']);
        print(gcf,'-dpdf',filename)
        clear filename
    end
end
close all

%% Compare partitions - adjusted mutual information

adjMI = double.empty;
ci = 0; % counting indices
for c=2:Nc % every column (g1)
for r=2:Nr % every row (g1)
    ci=ci+1;
    cons_labels{ci} = [cellData{r,1} ' ' cellData{1,c}];
    cii = 0;
    for cc=2:Nc % every column (g2)
    for rr=2:Nr % every row (g2)
        cii=cii+1;
        p1 = mrccPartition{r,c};           %    
        p2 = mrccPartition{rr,cc};         %  
        adjMI(ci,cii)=ami(p1,p2); 
        clear p1 p2
    end
    end
end
end
% zero negative ami
adjMI(adjMI<0)=0;

% compare partitions across groups
figure(...
        'units','inches',...
        'position',[1 1 7.5 6],...
        'paperpositionmode','auto');
%%------------------------%%
imagesc(adjMI); axis square
yticks(1:1:cii); yticklabels(cons_labels)
xticks(1:1:cii); xticklabels(cons_labels); xtickangle(45)
title(['Adjusted Mutual Information - ' grp])
colorbar
filename = fullfile(outdir, ['adj_mutual_info_' grp '_consensus_comparisons.pdf']);
print(gcf,'-dpdf',filename)  % SAVE FIGURE
clear filename

%% create matrix plots of consensus ordered covariance
for r=2:Nr % every row
    for c=2:Nc % every column
        [X,Y,idx_ord]=grid_communities(mrccPartition{r,c});
        figure(...
        'units','inches',...
        'position',[1 1 6 6],...
        'paperpositionmode','auto');
    imagesc(covMat{r,c}(idx_ord,idx_ord)); axis square
    colormap(cmap); caxis([-1 1])
    hold on
    plot(X,Y,'k','LineWidth',2)
    xticks(1:1:N); yticks(1:1:N);
    yticklabels(ROIlabels(idx_ord)); xticklabels(ROIlabels(idx_ord)); xtickangle(45)
    colorbar
    title('Community Ordered Pearson Covariance Matrix')
    sgtitle([grp ' ' cellData{r,1} ' ' cellData{1,c}])
    filename = fullfile(outdir, ['commOrd_covMat_' grp '_' cellData{r,1} '_' cellData{1,c} '.pdf']);
    print(gcf,'-dpdf',filename)
    clear filename
    end
end
close all

%% print a list of labels by community
cons_comm_roi = rNames;    cons_comm_roi(1,1:Nc) = cNames;    
for r = 2:Nr % loop over age
    for c = 2:Nc % loop over sex
        comm = mrccPartition{r,c}; % parition
        nc = unique(comm); % community indices
        cons_comm_roi{r,c} = cell.empty;
        for cm=1:length(nc) % for every index
            % extract roi names (rows) for each community (columns)
            tmp = ROIlabels(comm==nc(cm))';
            cons_comm_roi{r,c}(1:length(tmp),cm)= tmp;
            clear tmp
        end
    end
end

%% test differences across within levels (ie along rows/columns)

% pairwise testing of edge-level differences

% start with comparisons across columns
if Nc>2
    m = nchoosek(2:Nc,2);
    row_comparisons{1,1} = grp;
    for ii=2:Nr
    for jj=1:size(m,1)
        row_comparisons{1,2} = 'perm_ttest_pval';
        row_comparisons{end+1,1} = [rNames{ii} '_' cNames{m(jj,1)} ' vs. ' rNames{ii} '_' cNames{m(jj,2)}];  
        row_comparisons{end,2} = perm_ttest2_cov(cellData{ii,m(jj,1)},cellData{ii,m(jj,2)},10000);
    end
    end
else 
    row_comparisons = [];
end

% now comparisons across rows
if Nr>2
    m = nchoosek(2:Nr,2);
    col_comparisons{1,1} = grp;
    for ii=2:Nc
    for jj=1:size(m,1)
        col_comparisons{1,2} = 'perm_ttest_pval';
        col_comparisons{end+1,1} = [cNames{ii} '_' rNames{m(jj,1)} ' vs. ' cNames{ii} '_' rNames{m(jj,2)}];  
        col_comparisons{end,2} = perm_ttest2_cov(cellData{m(jj,1),ii},cellData{m(jj,2),ii},10000);
    end
    end
else
    col_comparisons = [];
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
np = length(pval);
% number of significant edges at permutation p-value threshold 
for p=1:np
if np==1
    plab = ['pval_' num2str(pval)];
else
    plab = ['pval_range_' num2str(max(pval)) '_' num2str(min(pval))];
end
cnt=0;
for i=2:size(covMat,2)
    for j=2:size(covMat,1)
        cnt=cnt+1;
        mask = logical(cov_permP{j,i}>pval(p));       % binary mask of NONsignificant edges
        mat = covMat{j,i};                                   % pull out covariance matrix
        mat(mask)=0;                                    % zero edged great than pval
        mat(logical(eye(length(mat))))=0;               % zero the diagonal (self connections)
        thr_cov(:,:,p,cnt) = mat;                       % store thresholded covariance matrix
        if p ==1 % store labels on first pass
            labels{cnt} = [rNames{j} '-' cNames{i}];% store labels
        end
        clear mask mat
    end
end

figure('Units','inches','Position',[1 1 8 4])
tiledlayout('flow')
for i = 1:cnt
    nexttile
    imagesc(thr_cov(:,:,p,i)); axis square
    caxis([-1 1]); colormap(cmap)
    xticks([]); yticks([])
    title(labels{i})
    colorbar
end
sgtitle([grp ' Thresholded Covariance at perm p < ' num2str(pval(p))])
filename = fullfile(outdir, ['thresh_covMats_' grp '_p' num2str(pval(p)) '.pdf']);
print(gcf,'-dpdf',filename)
clear filename

for i=1:cnt
    density(p,i) = density_und(thr_cov(:,:,p,i));
    [posStrength(:,p,i),negStrength(:,p,i),totalPosStrength(p,i),totalNegStrength(p,i)]...
        = strengths_und_sign(thr_cov(:,:,p,i));
    degree(:,p,i) = degrees_und(thr_cov(:,:,p,i));
    [~,cpsz] = get_components(double(thr_cov(:,:,p,i)~=0));
    max_component_size(p,i) = max(cpsz);
    min_component_size(p,i) = min(cpsz);
    [clust_coef(:,p,i)] = clustering_coef_wu_sign(thr_cov(:,:,p,i),3);
end

end

% set plotting type depending on whether its a single p-value or a range
if np>1
    pt = 0; else; pt = 1;
end

%% Density
figure
switch pt
    case 1
        bar(density); xticklabels(labels); ylabel('Network Density')
        title([grp ' group Density at P < ' num2str(pval)])
        xticks(1:1:length(density)); xticklabels(labels)
    case 0
        plot(density)
        legend(labels,'Location','best')
        xticks(1:1:np); xticklabels(pval); xlabel('permutation p-value threshold') 
        ylabel('Network Density')
        title([grp ' group Density across p-value range'])
end
filename = fullfile(outdir, ['netwDensity_' grp '_' plab '.pdf']);
print(gcf,'-dpdf',filename)
clear filename

%% Degree
switch pt
    case 1
        figure
        boxplot(squeeze(degree)); xticklabels(labels); ylabel('Degree Distribution')
        title([grp ' group Degree at P < ' num2str(pval)])
    case 0
        figure('Units','inches','Position',[1 1 8 4])
        idx=0;
        gvec = double.empty;
        degreevec = double.empty;
        for i=1:length(labels)
            d = degree(:,:,i);
            d = reshape(d,[],1);
            degreevec=vertcat(degreevec,d);
            clear d
            for j=1:np
                idx=idx+1;
                gv = zeros(size(degree,1),1)+idx;
                gvec = vertcat(gvec,gv);
                clear gv
            end
            idx=idx+1;
            m(1,i) = idx-(Nr-1);
            if i~=length(labels)
                degreevec(end+1,1)=nan;
                gvec(end+1,1) = idx;
            end
        end
        boxplot(degreevec,gvec)
        xticks(m); xticklabels(labels)
        ylabel('Degree Distribution')
        title({[grp ' group Degree across p-value range: '], num2str(pval)})
end
filename = fullfile(outdir, ['netwDegree_' grp '_' plab '.pdf']);
print(gcf,'-dpdf',filename)
clear filename

%% Strength
switch pt
    case 1
        figure
        subplot(2,2,1)
        bar(totalPosStrength); xticklabels(labels); ylabel('Sum Positive Nodal Strength')
        title('Positive Strength')
        subplot(2,2,2)
        bar(totalNegStrength); xticklabels(labels); ylabel('Sum Negative Nodal Strength')
        title('Negative Strength')
        subplot(2,2,3)
        boxplot(squeeze(posStrength)); xticklabels(labels); ylabel('Positive Nodal Strength Distribution')
        subplot(2,2,4)
        boxplot(squeeze(negStrength)); xticklabels(labels); ylabel('Negative Nodal Strength Distribution')
        sgtitle([grp ' group Nodal and Cumulative Strength at P < ' num2str(pval)])

        filename = fullfile(outdir, ['netwStrength_' grp '_' plab '.pdf']);
        print(gcf,'-dpdf',filename)
        clear filename
    case 0
        figure('Units','inches','Position',[1 1 8 4])
        subplot(1,2,1)
        plot(totalPosStrength)
        legend(labels,'Location','best')
        xticks(1:1:np); xticklabels(pval); xlabel('permutation p-value threshold') 
        ylabel('Sum Positive Strength')
        subplot(1,2,2)
        plot(totalNegStrength)
        legend(labels,'Location','best')
        xticks(1:1:np); xticklabels(pval); xlabel('permutation p-value threshold')
        ylabel('Sum Negative Strength')
        sgtitle({[grp ' group Cumulative Strength across p-value range:'],num2str(pval)})
        
        filename = fullfile(outdir, ['netwTotalStrength_' grp '_' plab '.pdf']);
        print(gcf,'-dpdf',filename)
        clear filename

        idx=0;
        psvec = double.empty;
        nsvec = double.empty;
        for i=1:length(labels)
            ps = posStrength(:,:,i); ps = reshape(ps,[],1);
            ns = negStrength(:,:,i); ns = reshape(ns,[],1);
            psvec=vertcat(psvec,ps); nsvec=vertcat(nsvec,ns);
            clear ps ns
            idx=idx+1+np;
            if i~=length(labels)
                psvec(end+1,1)=nan;
                nsvec(end+1,1)=nan;
            end
        end
        figure('Units','inches','Position',[1 1 8 8])
        subplot(2,2,1:2)
        boxplot(psvec,gvec)
        xticks(m); xticklabels(labels)
        ylabel('Positive Nodal Strength')
        subplot(2,2,3:4)
        boxplot(nsvec,gvec)
        xticks(m); xticklabels(labels)
        ylabel('Negative Nodal Strength')
        sgtitle({[grp ' group Strength across p-value range: '], num2str(pval)})

        filename = fullfile(outdir, ['netwNodalStrength_' grp '_' plab '.pdf']);
        print(gcf,'-dpdf',filename)
        clear filename
end
       
%% Largest Connectected Component
switch pt
    case 1
        figure
        subplot(1,2,1)
        bar(max_component_size); xticklabels(labels); ylabel('Number of Nodes')
        subplot(1,2,2)
        idx = find(min_component_size==1);
        if ~isempty(idx)
            hold on
            ylim([0 6]); yticks(1:1:6)
            xticks([])
            yticklabels(labels)
            for i=1:length(idx)
                text(0.1,idx(i),{'Contains fractured nodes of component size 1'})
            end
        end
        sgtitle([grp ' group Largest Connected Component at P < ' num2str(pval)])
    case 0
        figure('Units','inches','Position',[1 1 8 4])
        subplot(1,2,1)
        plot(max_component_size)
        legend(labels,'Location','best')
        xticks(1:1:np); xticklabels(pval); xlabel('permutation p-value threshold')
        ylabel('Number of nodes in largest component')
        subplot(1,2,2)
        idx = min_component_size==1;
        imagesc(idx)
        yticks(1:1:np); yticklabels(pval); ylabel('p-value threshold')
        xticks(1:1:length(labels)); xticklabels(labels)
        title('fractured node presence = 1')
        colorbar
        sgtitle([grp ' group Largest Connected Component across p-value range:'])
end
filename = fullfile(outdir, ['netwLargestComponent_' grp '_' plab '.pdf']);
print(gcf,'-dpdf',filename)
clear filename
%% Clustering Coefficient
switch pt
    case 1
        figure
        boxplot(squeeze(clust_coef)); xticklabels(labels); ylabel('Clusterring Coefficient')
        title([grp ' group Clustering Coefficient at P < ' num2str(pval)])
    case 0
        figure('Units','inches','Position',[1 1 8 4])
        idx=0;
      
        ccvec = double.empty;
        for i=1:length(labels)
            cc = clust_coef(:,:,i);
            cc = reshape(cc,[],1);
            ccvec=vertcat(ccvec,cc);
            clear cc
            idx=idx+1+np;
            if i~=length(labels)
                ccvec(end+1,1)=nan;
            end
        end
        boxplot(ccvec,gvec)
        ylim([0 1])
        xticks(m); xticklabels(labels)
        ylabel('Clustering Coefficient')
        title({[grp ' group Clustering Coefficient across p-value range: '], num2str(pval)})
end
filename = fullfile(outdir, ['netwClustCoef_' grp '_' plab '.pdf']);
print(gcf,'-dpdf',filename)
clear filename
%%
%-------------------------------------------------------------------------%
ver = length(dir(fullfile(pwd,['covariance_out_tier1_' grp '*'])));
if ver == 0
    fileout = ['covariance_out_tier1_' grp '.mat'];
else
    fileout = ['covariance_out_tier1_' grp '_run' num2str(ver+1) '.mat'];
end
if exist('covariates','var')
    save(fileout,'adjMI','agrMat','allPartitions',...
        'cellData','cNames','cons_comm_roi','cons_labels','cov_permP','cov_paramP',...
        'covMat','grp','mrccPartition','N','nc','Nc','Nr','R','rNames','ROIlabels','S',...
        'col_comparisons','row_comparisons',...
        'clust_coef','degree','density','labels','max_component_size','min_component_size',...
        'negStrength','posStrength','pval','thr_cov','totalPosStrength','totalNegStrength',...
        'covariates')
else
    save(fileout,'adjMI','agrMat','allPartitions',...
        'cellData','cNames','cons_comm_roi','cons_labels','cov_permP','cov_paramP',...
        'covMat','grp','mrccPartition','N','nc','Nc','Nr','R','rNames','ROIlabels','S',...
        'col_comparisons','row_comparisons',...
        'clust_coef','degree','density','labels','max_component_size','min_component_size',...
        'negStrength','posStrength','pval','thr_cov','totalPosStrength','totalNegStrength')
end
close all










