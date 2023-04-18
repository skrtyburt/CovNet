function covariance_analysis_tier2(g1t1_struct,g2t1_struct,cmask,pthresh)
narginchk(2,4)
% Input
%   g1t1_struct - output from the tier 1 script for group 1 loaded into a 
%                 matlab structure.
%   g2t1_struct - output from the tier 1 script for group 2 loaded into a
%                 matlab structure.
% optional:
%   cmask       - Binary matrix equal in size to data elements in the cellData
%                 that was used as tier1 input. It is used to isolate subsets of groups
%                 to be analyzed.
%   pval        - p-values for threholded networks examined in tier1 to be
%                 compared here.
%
% Author: Evgeny Chumin (2023), Indiana University
%%
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
%%
% find cells from which to compute summary values
lidx = find(cmask);

% set outout directory
%  CREATE NEW NAME FROM GROUP CONCATENATION
outroot = ['tier2_' g1t1_struct.grp '-' g2t1_struct.grp '_comparison'];

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
%%

%% pairwise comparison of edge distribution between g1 and g2 (UNTHRESHOLDED)
N = g1t1_struct.N;
Nr = g1t1_struct.Nr;
Nc = g1t1_struct.Nc;

idx = 0;
mask = logical(triu(ones(N,N),1));
for rw = 2:Nr
    for cl = 2:Nc
        idx=idx+1;
        grp{1,idx}=[g1t1_struct.covMat{rw,1} '-' g1t1_struct.covMat{1,cl}];
        tmp1 = g1t1_struct.covMat{rw,cl}(mask);
        tmp2 = g2t1_struct.covMat{rw,cl}(mask);
        rVec1(:,idx)= tmp1; rVec2(:,idx)= tmp2;
        clear tmp1 tmp2
    end
end
Nbins=round(size(rVec1,1)*.2);
comp=[g1t1_struct.grp ' vs. ' g2t1_struct.grp];
f = figure('units','inches','position',[1 1 8 6],'paperpositionmode','auto');
t = tiledlayout('flow');
for ix = 1:length(grp)
    nexttile
    histogram(rVec1(:,ix),Nbins,'DisplayStyle','bar','LineStyle','none','FaceColor','b','FaceAlpha',.1);hold on
    histogram(rVec2(:,ix),Nbins,'DisplayStyle','bar','LineStyle','none','FaceColor','r','FaceAlpha',.1);
    box off
    [~,p,ks]=kstest2(rVec1(:,ix),rVec2(:,ix));
    title([grp{ix} ': ks=' num2str(ks) ' p=' num2str(p)])
end
title(t,comp)
print(f,'-dpdf',[outroot '/' comp '_edge_histograms.pdf'])

%% pairwise comparisons of edge-level differences
idx = 0;
for rw = 2:Nr
    for cl = 2:Nc
        idx=idx+1;
        cdVec(idx,1) = g1t1_struct.cellData(rw,cl);
        cdVec(idx,2) = g2t1_struct.cellData(rw,cl);
    end
end
for ix=1:length(grp)
    [permp(:,:,ix), mcpermp(:,:,ix)] = perm_ttest2_cov(cdVec{ix,1},cdVec{ix,2},10000);
    rIX(:,ix)=logical(sum(permp(:,:,ix)<0.05));
    temp = permp(:,:,ix); temp=temp(mask);
    [pfdr]=FDR(temp,0.05);
    clear temp
end
rIX=logical(rIX);
f = figure('units','inches','position',[1 1 12 6],'paperpositionmode','auto');
t = tiledlayout('flow');
for ix = 1:length(grp)
    nexttile
    imagesc(permp(:,:,ix),[0 0.05]); axis square;
    colormap(flipud(hot)); colorbar
        hold on;
    spy(permp(:,:,ix)<pfdr,'o')
    spy(mcpermp(:,:,ix)<0.05,'c*')
    ylim([.5 N+.5]); xlim([.5 N+.5])
    tc=find(rIX(:,ix));
    ax=gca;
    xticks(tc); yticks(tc);
    xticklabels(g1t1_struct.ROIlabels(rIX(:,ix))); xtickangle(30)
    yticklabels(g1t1_struct.ROIlabels(rIX(:,ix)))
    ax.FontSize=6;
    title(grp{ix})
    clear tc
end
title(t,{comp, 'Edgewise difference testing. o-FDR<0.05 *-max statistic<0.05'})
print(f,'-dpdf',[outroot '/' comp '_edgewise_groupdiff.pdf'],'-fillpage')

%% compare nodal degree distributions
for p=1:length(pidx) % for every p-value
    pcmp = g1t1_struct.pval(pidx(p));
    Nbins=round(N/round(N*.1));
    mxidx = max([max(max(g1t1_struct.degree(:,pidx(p),:))) max(max(g2t1_struct.degree(:,pidx(p),:)))]);
    for ll=1:length(lidx)
        hmax(ll,1) = max(histcounts(g1t1_struct.degree(:,pidx(p),lidx(ll)),Nbins));
        hmax(ll,2) = max(histcounts(g2t1_struct.degree(:,pidx(p),lidx(ll)),Nbins));
    end
    hmax = max(max(hmax));
    f = figure('units','inches','position',[1 1 10.5 6],'PaperOrientation','landscape');
    t = tiledlayout('flow');
    for ll=1:length(lidx)
        nexttile
        dd1 = g1t1_struct.degree(:,pidx(p),lidx(ll));
        dd2 = g2t1_struct.degree(:,pidx(p),lidx(ll));
        
        histogram(dd1,'DisplayStyle','stairs','EdgeAlpha',.5,'LineWidth',2,'NumBins',Nbins)
        hold on
        histogram(dd2,'DisplayStyle','stairs','EdgeAlpha',.5,'LineWidth',2,'NumBins',Nbins)
        xlim([0 mxidx+1]); xticks(0:2:mxidx+1)
        ylim([0 hmax+1]); yticks(0:2:hmax+1)
        [~,ksp]=kstest2(dd1,dd2);
        xlabel('Degree'); ylabel('Frequency')
        title({grp{ll}, ['Degree Distributions: KS test p = ' num2str(ksp)]})
        legend({g1t1_struct.grp,g2t1_struct.grp},'Location','northeast')
        clear dd1 dd2
    end
    title(t,{[comp ' Degree Distributions'], ['network threshold p<' num2str(pcmp)]})
    filename = fullfile(outdir, ['Degree_dist_' comp '_netw_pthr' num2str(pcmp) '.pdf']);
    print(f,'-dpdf',filename)
    clear mxidx hmax filename
end
%% compare positive strength distributions
for p=1:length(pidx) % for every p-value
    mxidx = max([max(max(g1t1_struct.posStrength(:,pidx(p),:))) max(max(g2t1_struct.posStrength(:,pidx(p),:)))]);
    for ll=1:length(lidx)
        hmax(ll,1) = max(histcounts(g1t1_struct.posStrength(:,pidx(p),lidx(ll)),Nbins));
        hmax(ll,2) = max(histcounts(g2t1_struct.posStrength(:,pidx(p),lidx(ll)),Nbins));
    end
    hmax = max(max(hmax));
    f = figure('units','inches','position',[1 1 10.5 6],'PaperOrientation','landscape');
    t = tiledlayout('flow');
    for ll=1:length(lidx)
        nexttile
        dd1 = g1t1_struct.posStrength(:,pidx(p),lidx(ll));
        dd2 = g2t1_struct.posStrength(:,pidx(p),lidx(ll));

        histogram(dd1,'DisplayStyle','stairs','EdgeAlpha',.5,'LineWidth',2,'NumBins',Nbins)
        hold on
        histogram(dd2,'DisplayStyle','stairs','EdgeAlpha',.5,'LineWidth',2,'NumBins',Nbins)
        xlim([0 mxidx+1]); xticks(0:2:mxidx+1)
        ylim([0 hmax+1]); yticks(0:2:hmax+1)
        [~,ksp]=kstest2(dd1,dd2);
        xlabel('Positive Strength'); ylabel('Frequency')
        [tl,~]=title({grp{ll}, 'Positive Strength Distributions',['KS test p = ' num2str(ksp)]});
        tl.FontSize=10;
        legend({g1t1_struct.grp,g2t1_struct.grp},'Location','northeast')
        clear dd1 dd2
    end
    title(t,{[comp ' Positive Strength Distributions'], ['network threshold p<' num2str(pcmp)]})
    filename = fullfile(outdir, ['Pos_strength_' comp '_netw_pthr' num2str(pcmp) '.pdf']);
    % print(f,'-dpdf',filename)
    clear mxidx hmax filename 
end


%% compare negative strength distributions
    
    mxidx = max([max(max(g1t1_struct.negStrength(:,pidx(p),:))) max(max(g2t1_struct.negStrength(:,pidx(p),:)))]);
    for ll=1:length(lidx)
        hmax(ll,1) = max(histcounts(g1t1_struct.negStrength(:,pidx(p),lidx(ll)),nbins));
        hmax(ll,2) = max(histcounts(g2t1_struct.negStrength(:,pidx(p),lidx(ll)),nbins));
    end
    hmax = max(max(hmax));
    for ll=1:length(lidx)
        dd1 = g1t1_struct.negStrength(:,pidx(p),lidx(ll));
        lb1 = [g1t1_struct.grp '-' g1t1_struct.labels{lidx(ll)}];
        dd2 = g2t1_struct.negStrength(:,pidx(p),lidx(ll));
        lb2 = [g2t1_struct.grp '-' g2t1_struct.labels{lidx(ll)}];
        figure;
        histogram(dd1,'DisplayStyle','stairs','EdgeAlpha',.5,'LineWidth',2,'NumBins',nbins)
        hold on
        histogram(dd2,'DisplayStyle','stairs','EdgeAlpha',.5,'LineWidth',2,'NumBins',nbins)
        xlim([0 mxidx+1]); xticks(0:2:mxidx+1)
        ylim([0 hmax+1]); yticks(0:2:hmax+1)
        [~,ksp]=kstest2(dd1,dd2);
        xlabel('Negative Strength'); ylabel('Frequency')
        title([num2str(pcmp) ' thr: Negative Strength Distributions: KS test p = ' num2str(ksp)])
        legend({lb1,lb2},'Location','northeastoutside')

        filename = fullfile(outdir, ['neg_strength_dist_' lb1 '_v_' lb2 '_netw_pthr' num2str(pcmp) '.pdf']);
        print(gcf,'-dpdf',filename)
        clear filename dd1 dd2 lb1 lb2
    end
    clear mxidx hmax


%% compare clustering coefficient distributions

    
    mxidx = max([max(max(g1t1_struct.clust_coef(:,pidx(p),:))) max(max(g2t1_struct.clust_coef(:,pidx(p),:)))]);
    for ll=1:length(lidx)
        hmax(ll,1) = max(histcounts(g1t1_struct.clust_coef(:,pidx(p),lidx(ll)),nbins));
        hmax(ll,2) = max(histcounts(g2t1_struct.clust_coef(:,pidx(p),lidx(ll)),nbins));
    end
    hmax = max(max(hmax));
    for ll=1:length(lidx)
        dd1 = g1t1_struct.clust_coef(:,pidx(p),lidx(ll));
        lb1 = [g1t1_struct.grp '-' g1t1_struct.labels{lidx(ll)}];
        dd2 = g2t1_struct.clust_coef(:,pidx(p),lidx(ll));
        lb2 = [g2t1_struct.grp '-' g2t1_struct.labels{lidx(ll)}];
        figure;
        histogram(dd1,'DisplayStyle','stairs','EdgeAlpha',.5,'LineWidth',2,'NumBins',nbins)
        hold on
        histogram(dd2,'DisplayStyle','stairs','EdgeAlpha',.5,'LineWidth',2,'NumBins',nbins)
        xlim([0 mxidx+1]); xticks(0:2:mxidx+1)
        ylim([0 hmax+1]); yticks(0:2:hmax+1)
        [~,ksp]=kstest2(dd1,dd2);
        xlabel('Clustering Coefficient'); ylabel('Frequency')
        title([num2str(pcmp) ' thr: Clustering Coefficient Distributions: KS test p = ' num2str(ksp)])
        legend({lb1,lb2},'Location','northeastoutside')

        filename = fullfile(outdir, ['clust_coef_dist_' lb1 '_v_' lb2 '_netw_pthr' num2str(pcmp) '.pdf']);
        print(gcf,'-dpdf',filename)
        clear filename dd1 dd2 lb1 lb2
    end
    clear pcmp mxidx hmax
end