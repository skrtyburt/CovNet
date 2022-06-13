function covariance_analysis_tier2b(ci,label_ci,g1t1_struct,g2t1_struct,cmask,cmap,pthresh)

narginchk(4,7)
% check to make sure inputs are structures
if ~isstruct(g1t1_struct) || ~isstruct(g2t1_struct)
    fprintf(2,'input group tier1 outputs are not structures. Exiting...\n')
    return
end
% check to make sure row col names are the same among the two structures
if length(g1t1_struct.cNames)~=length(g2t1_struct.cNames) || ...
        length(g1t1_struct.rNames)~=length(g2t1_struct.rNames)
    fprintf(2,'Unequal number of rows/columns between group cells. Exiting...\n')
    return
end
% check node size vs size of community structure vector
if g1t1_struct.N ~= length(ci) || g2t1_struct.N ~= length(ci)
    fprintf(2,'Number of nodes does not match between inputs. Exiting...\n')
    return
end
N = g1t1_struct.N;
% make sure communities are linearly indexed
if sum(ci==1)==0 || length(unique(diff(unique(ci)))) > 1 
    % if index doesnt start with one or if there are gaps in index
    ci = fcn_relabel_partitions(ci);
end
% optional cmask should be the same size as number of DATA cells in the
% original input cell structure
if exist('cmask','var')
    [r,c]=size(cmask);
    if r ~= g1t1_struct.Nr-1 || c ~= g1t1_struct.Nc-1
        fprintf(2,'size of cmask does not match size data cells. Exiting...\n')
        return
    end
else
    cmask=ones(tier1_out_group1.Nr-1,tier1_out_group1.Nc-1);
end
%
if ~exist('cmap','var')
    cmap = turbo(64);
end
if length(g1t1_struct.pval)~=length(g2t1_struct.pval)
    fprintf(2,'Unequal number of p-value thresholds between groups.\n')
    fprintf(2,'Make sure both groups were ran with the same p-value inputs. Exiting...\n')
    return
elseif sum(ismember(g1t1_struct.pval,g2t1_struct.pval))~=length(g1t1_struct.pval)
    fprintf(2,'P-value thresholds do not match between groups. Exiting...\n')
    return
end
if exist('pthresh','var')
    pidx = find(ismember(g1t1_struct.pval,pthresh));
else
    pidx = 1:1:length(g1t1_struct.pval);
end
if isempty(pidx)
    fprintf(2,'Input pthresh values not found in group data pvalue thresholds. Exiting...\n')
    return
end

% find cells from which to compute summary values
lidx = find(cmask);
cmask = vertcat(zeros(1,g1t1_struct.Nc), horzcat(zeros(g1t1_struct.Nr-1,1),cmask));
[rw,cl]=find(cmask);
% find unique community labels
ciu = unique(ci);
% set outout directory
outdir = ['../' label_ci 'ord_mod_comparison'];
if ~exist(outdir,'dir')
    mkdir(outdir)
end

% plot and save the reference matrix ordered by the reference partition
% use the respective group title and ref partition
%% create matrix plots of consensus ordered covariance
[X,Y,idx_ord]=grid_communities(ci);
for r=1:length(rw) % every row
    figure(...
    'units','inches',...
    'position',[1 1 8.5 4],...
    'paperpositionmode','auto');
    subplot(1,2,1)
    imagesc(g1t1_struct.covMat{rw(r),cl(r)}(idx_ord,idx_ord)); axis square
    colormap(cmap); caxis([-1 1])
    hold on
    plot(X,Y,'k','LineWidth',2)
    xticks(1:1:N); yticks(1:1:N);
    yticklabels(g1t1_struct.ROIlabels(idx_ord)); xticklabels(g1t1_struct.ROIlabels(idx_ord)); xtickangle(45)
    colorbar
    title([g1t1_struct.grp ' ' g1t1_struct.cellData{rw(r),1} ' ' g1t1_struct.cellData{1,cl(r)}])

    subplot(1,2,2)
    imagesc(g2t1_struct.covMat{rw(r),cl(r)}(idx_ord,idx_ord)); axis square
    colormap(cmap); caxis([-1 1])
    hold on
    plot(X,Y,'k','LineWidth',2)
    xticks(1:1:N); yticks(1:1:N);
    yticklabels(g2t1_struct.ROIlabels(idx_ord)); xticklabels(g2t1_struct.ROIlabels(idx_ord)); xtickangle(45)
    colorbar
    title([g2t1_struct.grp ' ' g2t1_struct.cellData{rw(r),1} ' ' g2t1_struct.cellData{1,cl(r)}])
    
    sgtitle([label_ci ' ordered covariance matrices (Unthresholded)'],'Interpreter','none')

    filename = fullfile(outdir, [label_ci '_Ord_covMat_' g1t1_struct.cellData{rw(r),1} '_' g1t1_struct.cellData{1,cl(r)} '.pdf']);
    print(gcf,'-dpdf',filename)
    clear filename
end
close all

%% compare nodal degree distributions


for p=1:length(pidx)
    % loop over p-value thresholds
    pcmp = g1t1_struct.pval(pidx(p));
    for cmi = 1:length(ciu)
        % loop over communities
        for ll=1:length(lidx)
            % loop over groups in cmask
            g1m
            g2m






for p=1:length(pidx)
    pcmp = g1t1_struct.pval(pidx(p));
    for ll=1:length(lidx)
        dd1 = g1t1_struct.degree(:,pidx(p),lidx(ll));
        lb1 = [g1t1_struct.grp ' ' g1t1_struct.labels{lidx(ll)}];
        dd2 = g2t1_struct.degree(:,pidx(p),lidx(ll));
        lb2 = [g2t1_struct.grp ' ' g2t1_struct.labels{lidx(ll)}];
        figure;
        histogram(dd1,'DisplayStyle','stairs','EdgeAlpha',.5)
        hold on
        histogram(dd2,'DisplayStyle','stairs','EdgeAlpha',.5)
        [~,ksp]=kstest2(dd1,dd2);
        xlabel('Degree'); ylabel('Frequency')
        title(['Degree Distributions: KS test p = ' num2str(ksp)])
        legend({lb1,lb2},'Location','best')

        filename = fullfile(outdir, ['Degree_dist_' lb1 '_v_' lb2 'netw_pthr' num2str(pcmp) '.pdf']);
        print(gcf,'-dpdf',filename)
        clear filename dd1 dd2 lb1 lb2
    end
    clear pcmp
end
    
%compare within module degree

%% compare positive strength distributions
for p=1:length(pidx)
    pcmp = g1t1_struct.pval(pidx(p));
    for ll=1:length(lidx)
        dd1 = g1t1_struct.posStrength(:,pidx(p),lidx(ll));
        lb1 = [g1t1_struct.grp '-' g1t1_struct.labels{lidx(ll)}];
        dd2 = g2t1_struct.posStrength(:,pidx(p),lidx(ll));
        lb2 = [g2t1_struct.grp '-' g2t1_struct.labels{lidx(ll)}];
        figure;
        histogram(dd1,'DisplayStyle','stairs','EdgeAlpha',.5)
        hold on
        histogram(dd2,'DisplayStyle','stairs','EdgeAlpha',.5)
        [~,ksp]=kstest2(dd1,dd2);
        xlabel('Positive Strength'); ylabel('Frequency')
        title(['Positive Strength Distributions: KS test p = ' num2str(ksp)])
        legend({lb1,lb2},'Location','best')

        filename = fullfile(outdir, ['pos_strength_dist_' lb1 '_v_' lb2 'netw_pthr' num2str(pcmp) '.pdf']);
        print(gcf,'-dpdf',filename)
        clear filename dd1 dd2 lb1 lb2
    end
    clear pcmp
end
%% compare negative strength distributions
for p=1:length(pidx)
    pcmp = g1t1_struct.pval(pidx(p));
    for ll=1:length(lidx)
        dd1 = g1t1_struct.negStrength(:,pidx(p),lidx(ll));
        lb1 = [g1t1_struct.grp '-' g1t1_struct.labels{lidx(ll)}];
        dd2 = g2t1_struct.negStrength(:,pidx(p),lidx(ll));
        lb2 = [g2t1_struct.grp '-' g2t1_struct.labels{lidx(ll)}];
        figure;
        histogram(dd1,'DisplayStyle','stairs','EdgeAlpha',.5)
        hold on
        histogram(dd2,'DisplayStyle','stairs','EdgeAlpha',.5)
        [~,ksp]=kstest2(dd1,dd2);
        xlabel('Negative Strength'); ylabel('Frequency')
        title(['Negative Strength Distributions: KS test p = ' num2str(ksp)])
        legend({lb1,lb2},'Location','best')

        filename = fullfile(outdir, ['neg_strength_dist_' lb1 '_v_' lb2 'netw_pthr' num2str(pcmp) '.pdf']);
        print(gcf,'-dpdf',filename)
        clear filename dd1 dd2 lb1 lb2
    end
    clear pcmp
end

%% compare clustering coefficient distributions
for p=1:length(pidx)
    pcmp = g1t1_struct.pval(pidx(p));
    for ll=1:length(lidx)
        dd1 = g1t1_struct.clust_coef(:,pidx(p),lidx(ll));
        lb1 = [g1t1_struct.grp '-' g1t1_struct.labels{lidx(ll)}];
        dd2 = g2t1_struct.clust_coef(:,pidx(p),lidx(ll));
        lb2 = [g2t1_struct.grp '-' g2t1_struct.labels{lidx(ll)}];
        figure;
        histogram(dd1,'DisplayStyle','stairs','EdgeAlpha',.5)
        hold on
        histogram(dd2,'DisplayStyle','stairs','EdgeAlpha',.5)
        [~,ksp]=kstest2(dd1,dd2);
        xlabel('Clustering Coefficient'); ylabel('Frequency')
        title(['Clustering Coefficient Distributions: KS test p = ' num2str(ksp)])
        legend({lb1,lb2},'Location','best')

        filename = fullfile(outdir, ['clust_coef_dist_' lb1 '_v_' lb2 'netw_pthr' num2str(pcmp) '.pdf']);
        print(gcf,'-dpdf',filename)
        clear filename dd1 dd2 lb1 lb2
    end
    clear pcmp
end


%[Mden]=module_density_new(tier1_out_group1.covMat{rw(i),cl(i)},ci,1,'median');








