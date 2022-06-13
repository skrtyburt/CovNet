function covariance_analysis_tier2a(g1t1_struct,g2t1_struct,cmask,pthresh)

% Group comparisons of tier1 outputs within factor levels.

narginchk(2,4)
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
% check number of nodes
if g1t1_struct.N ~= g2t1_struct.N
    fprintf(2,'Number of nodes does not match between inputs. Exiting...\n')
    return
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

% set outout directory
%  CREATE NEW NAME FROM GROUP CONCATENATION
outroot = [g1t1_struct.grp '-' g2t1_struct.grp '_comparison'];

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

%% compare nodal degree distributions
for p=1:length(pidx)
    pcmp = g1t1_struct.pval(pidx(p));
    for ll=1:length(lidx)
        dd1 = g1t1_struct.degree(:,pidx(p),lidx(ll));
        lb1 = [g1t1_struct.grp ' ' g1t1_struct.labels{lidx(ll)}];
        dd2 = g2t1_struct.degree(:,pidx(p),lidx(ll));
        lb2 = [g2t1_struct.grp ' ' g2t1_struct.labels{lidx(ll)}];
        figure;
        histogram(dd1,'DisplayStyle','stairs','EdgeAlpha',.5,'LineWidth',2,'NumBins',5)
        hold on
        histogram(dd2,'DisplayStyle','stairs','EdgeAlpha',.5,'LineWidth',2,'NumBins',5)
        [~,ksp]=kstest2(dd1,dd2);
        xlabel('Degree'); ylabel('Frequency')
        title([num2str(pcmp) ' thr: Degree Distributions: KS test p = ' num2str(ksp)])
        legend({lb1,lb2},'Location','best')

        filename = fullfile(outdir, ['Degree_dist_' lb1 '_v_' lb2 '_netw_pthr' num2str(pcmp) '.pdf']);
        print(gcf,'-dpdf',filename)
        clear filename dd1 dd2 lb1 lb2
    end
    clear pcmp
end

%% compare positive strength distributions
for p=1:length(pidx)
    pcmp = g1t1_struct.pval(pidx(p));
    for ll=1:length(lidx)
        dd1 = g1t1_struct.posStrength(:,pidx(p),lidx(ll));
        lb1 = [g1t1_struct.grp '-' g1t1_struct.labels{lidx(ll)}];
        dd2 = g2t1_struct.posStrength(:,pidx(p),lidx(ll));
        lb2 = [g2t1_struct.grp '-' g2t1_struct.labels{lidx(ll)}];
        figure;
        histogram(dd1,'DisplayStyle','stairs','EdgeAlpha',.5,'LineWidth',2,'NumBins',5)
        hold on
        histogram(dd2,'DisplayStyle','stairs','EdgeAlpha',.5,'LineWidth',2,'NumBins',5)
        [~,ksp]=kstest2(dd1,dd2);
        xlabel('Positive Strength'); ylabel('Frequency')
        title([num2str(pcmp) ' thr: Positive Strength Distributions: KS test p = ' num2str(ksp)])
        legend({lb1,lb2},'Location','best')

        filename = fullfile(outdir, ['pos_strength_dist_' lb1 '_v_' lb2 '_netw_pthr' num2str(pcmp) '.pdf']);
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
        histogram(dd1,'DisplayStyle','stairs','EdgeAlpha',.5,'LineWidth',2,'NumBins',5)
        hold on
        histogram(dd2,'DisplayStyle','stairs','EdgeAlpha',.5,'LineWidth',2,'NumBins',5)
        [~,ksp]=kstest2(dd1,dd2);
        xlabel('Negative Strength'); ylabel('Frequency')
        title([num2str(pcmp) ' thr: Negative Strength Distributions: KS test p = ' num2str(ksp)])
        legend({lb1,lb2},'Location','best')

        filename = fullfile(outdir, ['neg_strength_dist_' lb1 '_v_' lb2 '_netw_pthr' num2str(pcmp) '.pdf']);
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
        histogram(dd1,'DisplayStyle','stairs','EdgeAlpha',.5,'LineWidth',2,'NumBins',5)
        hold on
        histogram(dd2,'DisplayStyle','stairs','EdgeAlpha',.5,'LineWidth',2,'NumBins',5)
        [~,ksp]=kstest2(dd1,dd2);
        xlabel('Clustering Coefficient'); ylabel('Frequency')
        title([num2str(pcmp) ' thr: Clustering Coefficient Distributions: KS test p = ' num2str(ksp)])
        legend({lb1,lb2},'Location','best')

        filename = fullfile(outdir, ['clust_coef_dist_' lb1 '_v_' lb2 '_netw_pthr' num2str(pcmp) '.pdf']);
        print(gcf,'-dpdf',filename)
        clear filename dd1 dd2 lb1 lb2
    end
    clear pcmp
end
close all
end
