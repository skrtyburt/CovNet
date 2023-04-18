function covariance_analysis_tier1(cellData,ROIlabels,cmap,pval)
narginchk(4,4)
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

% append repeated runs numerically
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
        [covMat{r,c}, cov_paramP{r,c}]=covariance(cellData{r,c},1,ROIlabels,[-1 1],cmap);
        sgtitle([grp ' ' cellData{r,1} ' ' cellData{1,c}])
        filename = fullfile(outdir, ['covMat_' grp '_' cellData{r,1} '_' cellData{1,c} '.pdf']);
        print(gcf,'-dpdf',filename)
        clear filename
    end
end
close all

% Parametric p-values of correlation significance are saved for
% completeness, but the recommended method for identifying signifant edge
% covariance is permutation testing.
%% Permutation testing of significant edge covariance 
% 10,000 permutations
cov_permP = rNames;    cov_permP(1,1:Nc) = cNames;    % permutation p value matrices
for r=2:Nr % every row
    for c=2:Nc % every column
        [~,~,cov_permP{r,c}] = randshiftnull_cov(cellData{r,c},10000);
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
colormap(cmap(33:end,:))
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

%% test differences across within leves (ie along rows/columns)
disp('Pairwise Permutation Testing of Covariance Matrices')
PairwisePermTests.alongrows = cell.empty;
nc = Nc-1;
ncomp = nchoosek(1:1:nc,2);
counter=0;
for r=2:Nr
    for cmp=1:size(ncomp,1)
        counter=counter+1;
        PairwisePermTests.alongrows{counter,1} = [rNames{r} '_' cNames{ncomp(cmp,1)+1} ' vs. ' rNames{r} '_' cNames{ncomp(cmp,2)+1}];
        PairwisePermTests.alongrows{counter,2} = perm_ttest2_cov(cellData{r,ncomp(cmp,1)+1},cellData{r,ncomp(cmp,2)+1},10000);
    end
end
clear counter
PairwisePermTests.alongcols = cell.empty;
nr = Nr-1;
ncomp = nchoosek(1:1:nr,2);
counter=0;
for c=2:Nc
    for cmp=1:size(ncomp,1)
        counter=counter+1;
        PairwisePermTests.alongcols{counter,1} = [cNames{c} '_' rNames{ncomp(cmp,1)+1} ' vs. ' cNames{c} '_' rNames{ncomp(cmp,2)+1}];
        PairwisePermTests.alongcols{counter,2} = perm_ttest2_cov(cellData{ncomp(cmp,1)+1,c},cellData{ncomp(cmp,2)+1,c},10000);
    end
end
clear counter
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for every selected p-value threshold
np = length(pval);
%% Create the thresholded covaraince matrices
for p=1:np
    cnt=0;
    for r=2:Nr
        for c=2:Nc
            cnt=cnt+1;
            mask = logical(cov_permP{r,c}>pval(p));     % binary mask of NONsignificant edges
            mat = covMat{r,c};                          % pull out covariance matrix
            if p==1
                cov(:,:,cnt) = mat;                     % save covmat for later analysis
            end
            mat(mask)=0;                                % zero edges greater than pval
            mat(logical(eye(length(mat))))=0;           % zero the diagonal (self connections)
            thr_cov(:,:,p,cnt) = mat;                   % store thresholded covariance matrix
            if p ==1 % store labels on first pass
                labels{cnt} = [rNames{r} '-' cNames{c}];% store labels
            end
            clear mask mat
        end
    end
%% Plot the thresholded covariance matrices
    figure('Units','inches','Position',[1 1 8 4])
    tiledlayout('flow')
    for i = 1:cnt
        nexttile
        imagesc(thr_cov(:,:,p,i)); axis square
        caxis([-1 1]); colormap(cmap)
        xticks(1:1:N); yticks(1:1:N);
        yticklabels(ROIlabels); xticklabels(ROIlabels); xtickangle(45)
        title(labels{i})
        colorbar
    end
    sgtitle({[grp ' Thresholded Covariance at perm p < ' num2str(pval(p))], ' '})
    filename = fullfile(outdir, ['thresh_covMats_' grp '_p' num2str(pval(p)) '.pdf']);
    print(gcf,'-dpdf',filename)
    clear filename
end

%% Compute global and regional network metrics
[~,~,np,ng]=size(thr_cov);
% build p-value labels
for p=1:np
    plabs{p}=['p<' num2str(pval(p))];
end
p_range = [num2str(max(pval)) '-' num2str(min(pval))];
%% Density
% only applies to thresholded networks
for p=1:np
    for g=1:ng
        density(p,g) = density_und(thr_cov(:,:,p,g));
    end
end
f=figure(...
        'units','inches',...
        'position',[1 1 5 4],...
        'paperpositionmode','auto');
bar(density); xticklabels(plabs); ylim([0 1]); ylabel('Network Density')
legend(labels,'Location','eastoutside')
filename = fullfile(outdir, ['netwDensity_' grp '_' p_range '.pdf']);
print(f,'-dpdf',filename)
clear filename f

%% Degree
% only applies to thresholded networks
for p=1:np
    for g=1:ng
        degree(:,p,g) = degrees_und(thr_cov(:,:,p,g));
    end
end

f = figure('Units','inches','Position',[1 1 8 4]);
degvec=double.empty;
deggrp=double.empty;
idx=1;
boxlabs=cell.empty;
for p=1:np
    for g=1:length(labels)
        degvec = vertcat(degvec,degree(:,p,g));
        deggrp = vertcat(deggrp,zeros(N,1)+idx);
        idx=idx+1;
    end
    degvec(end+1,1)=nan;
    deggrp(end+1,1)=idx;
    idx=idx+1;
    boxlabs = horzcat(boxlabs,labels,' ');
end

boxplot(degvec,deggrp)
xticklabels(boxlabs)
ylabel('Degree Distribution'); ylim([0 N])
box off
for p=1:np
    text((length(labels)+1)/2+length(labels)*(p-1),N-1,plabs{p})
end
title({[grp ' group Degree across p-value range: '], p_range})

filename = fullfile(outdir, ['netwDegree_' grp '_' p_range '.pdf']);
print(f,'-dpdf',filename)
clear filename f

%% Strength
% applies to both unthresholded and thresholded networks
% compute from unthresholded networks
for g=1:ng
    [posStrength(:,1,g),negStrength(:,1,g),totalPosStrength(1,g),totalNegStrength(1,g)]...
        = strengths_und_sign(cov(:,:,g));
end
% Compute from thresholded networks
for p=1:np
    for g=1:ng
        [posStrength(:,p+1,g),negStrength(:,p+1,g),totalPosStrength(p+1,g),totalNegStrength(p+1,g)]...
            = strengths_und_sign(thr_cov(:,:,p,g));
    end
end

% Format for boxplots
posSvec=double.empty;
negSvec=double.empty;
STRgrp=double.empty;
idx=1;
boxSTRlabs=cell.empty;
for p=1:np+1
    for g=1:length(labels)
        posSvec = vertcat(posSvec,posStrength(:,p,g));
        negSvec = vertcat(negSvec,negStrength(:,p,g));
        STRgrp = vertcat(STRgrp,zeros(N,1)+idx);
        idx=idx+1;
    end
    posSvec(end+1,1)=nan;
    negSvec(end+1,1)=nan;
    STRgrp(end+1,1)=idx;
    idx=idx+1;
    boxSTRlabs = horzcat(boxSTRlabs,labels,' ');
end
STRplabs=horzcat({'unthesholded'},plabs);

% PLOT: Nodal positive strength
f = figure('Units','inches','Position',[1 1 8 4]);
upr = max(posSvec)+.1*max(posSvec);
boxplot(posSvec,STRgrp)
xticklabels(boxSTRlabs)
ylabel('Nodal Postive Strengths'); ylim([0 upr])
box off
for p=1:np+1
    text((length(labels)+1)/2+length(labels)*(p-1),upr-(.02*upr),STRplabs{p})
end
title({[grp ' Nodal Positive Strength across p-value range: '], p_range})
filename = fullfile(outdir, ['nodePosStrength_' grp '_' p_range '.pdf']);
print(f,'-dpdf',filename)
clear filename f

% PLOT: Nodal negative strength
f = figure('Units','inches','Position',[1 1 8 4]);
upr = max(negSvec)+.1*max(negSvec);
boxplot(negSvec,STRgrp)
xticklabels(boxSTRlabs)
ylabel('Nodal Negative Strengths'); ylim([0 upr])
box off
for p=1:np+1
    text((length(labels)+1)/2+length(labels)*(p-1),upr-(.02*upr),STRplabs{p})
end
title({[grp ' Nodal Negative Strength across p-value range: '], p_range})
filename = fullfile(outdir, ['nodeNegStrength_' grp '_' p_range '.pdf']);
print(f,'-dpdf',filename)
clear filename f

% PLOT: Total positive Strength
f=figure(...
        'units','inches',...
        'position',[1 1 5 4],...
        'paperpositionmode','auto');
bar(totalPosStrength); xticklabels(STRplabs); ylabel('Total Positive Strength')
legend(labels,'Location','eastoutside')
filename = fullfile(outdir, ['totPosStrength_' grp '_' p_range '.pdf']);
print(f,'-dpdf',filename)
clear filename f

% PLOT: Total negative Strength
f=figure(...
        'units','inches',...
        'position',[1 1 5 4],...
        'paperpositionmode','auto');
bar(totalNegStrength); xticklabels(STRplabs); ylabel('Total Negative Strength')
legend(labels,'Location','eastoutside')
filename = fullfile(outdir, ['totNegStrength_' grp '_' p_range '.pdf']);
print(f,'-dpdf',filename)
clear filename f
       
%% Largest Connectected Component
% applies only to thresholded networks to find disconnected nodes/components
for p=1:np
    for g=1:ng
        [~,cpsz] = get_components(double(thr_cov(:,:,p,g)~=0));
        max_component_size(p,g) = max(cpsz);
        min_component_size(p,g) = min(cpsz);
        clear cpsz
    end
end
f=figure(...
        'units','inches',...
        'position',[1 1 5 4],...
        'paperpositionmode','auto');
bar(min_component_size); xticklabels(plabs); ylim([0 N]); ylabel('Size of Smallest Connected Component')
legend(labels,'Location','eastoutside')
title({['Y < ' num2str(N) ' = Fractured Network'],...
    'Y ==1 = Disconnected Node'})

filename = fullfile(outdir, ['minConnectedComponent_' grp '_' p_range '.pdf']);
print(f,'-dpdf',filename)
clear filename f

%% Clustering Coefficient
% The clustering coefficient variant used here can be applied to fully
% connected and thresholded network. 
for g=1:ng
    [clust_coef(:,1,g)] = clustering_coef_wu_sign(cov(:,:,g),3);
end
% Compute from thresholded networks
for p=1:np
    for g=1:ng
        [clust_coef(:,p+1,g)] = clustering_coef_wu_sign(thr_cov(:,:,p,g),3);
    end
end

% Format for boxplots
clustvec=double.empty;
idx=1;
for p=1:np+1
    for g=1:length(labels)
        clustvec = vertcat(clustvec,clust_coef(:,p,g));
        idx=idx+1;
    end
    clustvec(end+1,1)=nan;
    idx=idx+1;
end

% PLOT: clustering coefficient
f = figure('Units','inches','Position',[1 1 8 4]);
boxplot(clustvec,STRgrp)
xticklabels(boxSTRlabs)
ylabel('Clustering Coefficient'); ylim([0 1])
box off
for p=1:np+1
    text((length(labels)+1)/2+length(labels)*(p-1),.97,STRplabs{p})
end
title({[grp ' Clustering Coefficient across p-value range: '], p_range})
filename = fullfile(outdir, ['ClustCoeff_' grp '_' p_range '.pdf']);
print(f,'-dpdf',filename)
clear filename f

%%
%-------------------------------------------------------------------------%
ver = length(dir(fullfile(pwd,['covariance_out_tier1_' grp '*'])));
if ver == 0
    fileout = ['covariance_out_tier1_' grp '.mat'];
else
    fileout = ['covariance_out_tier1_' grp '_run' num2str(ver+1) '.mat'];
end
save(fileout,'adjMI','agrMat','allPartitions',...
    'cellData','cNames','cons_comm_roi','cons_labels','cov_permP','cov_paramP',...
    'covMat','grp','mrccPartition','N','nc','Nc','Nr','R','rNames','ROIlabels','S',...
    'PairwisePermTests',...
    'clust_coef','degree','density','labels','max_component_size','min_component_size',...
    'negStrength','posStrength','pval','thr_cov','totalPosStrength','totalNegStrength')
close all










