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
labels = g1t1_struct.labels(lidx);
% find unique community labels
ciu = unique(ci);
% set outout directory
outroot = ['tier2b_' label_ci 'ord_mod_comparison'];

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

%% Concatenate and plot group joing AMI matrix

adjMI = zeros(length(rw)*2);
cidx=0;
for r = 1:length(rw)
        cidx = cidx+1;
        gci(:,cidx) = g1t1_struct.mrccPartition{rw(r),cl(r)};
        gci(:,cidx+length(rw)) = g2t1_struct.mrccPartition{rw(r),cl(r)};
        
        gci_label{1,cidx} = [g1t1_struct.grp '-' g1t1_struct.rNames{rw(r)} '-' g1t1_struct.cNames{cl(r)}];
        gci_label{1,cidx+length(rw)} = [g2t1_struct.grp '-' g2t1_struct.rNames{rw(r)} '-' g2t1_struct.cNames{cl(r)}];
end
for r2 = 1:size(gci,2)-1
    for r3 = r2+1:size(gci,2)
        adjMI(r2,r3) = ami(gci(:,r2),gci(:,r3));
    end
end
adjMI = adjMI + adjMI';
adjMI(adjMI<0)=0;
diag = logical(eye(length(adjMI)));
adjMI(diag)=0;

% compare partitions across groups
figure(...
        'units','inches',...
        'position',[1 1 7.5 6],...
        'paperpositionmode','auto');
%%------------------------%%
imagesc(adjMI); axis square
yticks(1:1:length(gci_label)); yticklabels(gci_label)
xticks(1:1:length(gci_label)); xticklabels(gci_label); xtickangle(45)
title('Adjusted Mutual Information')
colorbar; colormap(cmap(length(cmap)/2:end,:))
filename = fullfile(outdir, 'adj_mutual_info_comparisons.pdf');
print(gcf,'-dpdf',filename)  % SAVE FIGURE
clear filename

%% compare within module degree

for p=1:length(pidx)
    % loop over p-value thresholds
    pcmp = g1t1_struct.pval(pidx(p));
    for cmi = 1:length(ciu)
        cidx = ci==ciu(cmi);
        % loop over communities
        for ll=1:length(lidx)
            % loop over groups in cmask
            g1_deg(:,ll) = g1t1_struct.degree(cidx,pidx(p),lidx(ll));
            g2_deg(:,ll) = g2t1_struct.degree(cidx,pidx(p),lidx(ll));
        end
        
        figure(...
        'units','inches',...
        'position',[1 1 7.5 4],...
        'paperpositionmode','auto');
    
        tiledlayout('flow')
        nexttile
        boxplot(g1_deg,labels)
        ylabel('Nodal Degree')
        title({' ',g1t1_struct.grp})
        
        nexttile
        boxplot(g2_deg,labels)
        title({' ',g2t1_struct.grp})
        
        sgtitle({[label_ci ' ref-partition w/in module degree'],['network thresh p<' num2str(pcmp) ' module ' num2str(ciu(cmi))]},'Interpreter','none')
        
        clear g1_deg g2_deg
        
        filename = fullfile(outdir, [label_ci '_refPart_within_degree_' g1t1_struct.grp '_vs_' g2t1_struct.grp '_thr' num2str(pcmp) '_m' num2str(ciu(cmi)) '.pdf']);
        print(gcf,'-dpdf',filename)
        clear filename
    end
    clear pcmp
    
end

%% compare within module positive strength

for p=1:length(pidx)
    % loop over p-value thresholds
    pcmp = g1t1_struct.pval(pidx(p));
    for cmi = 1:length(ciu)
        cidx = ci==ciu(cmi);
        % loop over communities
        for ll=1:length(lidx)
            % loop over groups in cmask
            g1_pstr(:,ll) = g1t1_struct.posStrength(cidx,pidx(p),lidx(ll));
            g2_pstr(:,ll) = g2t1_struct.posStrength(cidx,pidx(p),lidx(ll));
        end
        
        figure(...
        'units','inches',...
        'position',[1 1 7.5 4],...
        'paperpositionmode','auto');
    
        tiledlayout('flow')
        nexttile
        boxplot(g1_pstr,labels)
        ylabel('Nodal Positive Strength')
        title({' ',g1t1_struct.grp})
        
        nexttile
        boxplot(g2_pstr,labels)
        title({' ',g2t1_struct.grp})
        
        sgtitle({[label_ci ' ref-partition w/in module nodal pos strength'],['network thresh p<' num2str(pcmp) ' module ' num2str(ciu(cmi))]},'Interpreter','none')
        
        clear g1_pstr g2_pstr
        
        filename = fullfile(outdir, [label_ci '_refPart_within_posStrength_' g1t1_struct.grp '_vs_' g2t1_struct.grp '_thr' num2str(pcmp) '_m' num2str(ciu(cmi)) '.pdf']);
        print(gcf,'-dpdf',filename)
        clear filename
    end
    clear pcmp
    
end

%% compare within module negative strength

for p=1:length(pidx)
    % loop over p-value thresholds
    pcmp = g1t1_struct.pval(pidx(p));
    for cmi = 1:length(ciu)
        cidx = ci==ciu(cmi);
        % loop over communities
        for ll=1:length(lidx)
            % loop over groups in cmask
            g1_nstr(:,ll) = g1t1_struct.negStrength(cidx,pidx(p),lidx(ll));
            g2_nstr(:,ll) = g2t1_struct.negStrength(cidx,pidx(p),lidx(ll));
        end
        
        figure(...
        'units','inches',...
        'position',[1 1 7.5 4],...
        'paperpositionmode','auto');
    
        tiledlayout('flow')
        nexttile
        boxplot(g1_nstr,labels)
        ylabel('Nodal Negative Strength')
        title({' ',g1t1_struct.grp})
        
        nexttile
        boxplot(g2_nstr,labels)
        title({' ',g2t1_struct.grp})
        
        sgtitle({[label_ci ' ref-partition w/in module nodal neg strength'],['network thresh p<' num2str(pcmp) ' module ' num2str(ciu(cmi))]},'Interpreter','none')
        
        clear g1_nstr g2_nstr
        
        filename = fullfile(outdir, [label_ci '_refPart_within_negStrength_' g1t1_struct.grp '_vs_' g2t1_struct.grp '_thr' num2str(pcmp) '_m' num2str(ciu(cmi)) '.pdf']);
        print(gcf,'-dpdf',filename)
        clear filename
        
    end
    clear pcmp
    
end

%% compare within module clustering coefficient

for p=1:length(pidx)
    % loop over p-value thresholds
    pcmp = g1t1_struct.pval(pidx(p));
    for cmi = 1:length(ciu)
        cidx = ci==ciu(cmi);
        % loop over communities
        for ll=1:length(lidx)
            % loop over groups in cmask
            g1_cc(:,ll) = g1t1_struct.clust_coef(cidx,pidx(p),lidx(ll));
            g2_cc(:,ll) = g2t1_struct.clust_coef(cidx,pidx(p),lidx(ll));
        end
        
        figure(...
        'units','inches',...
        'position',[1 1 7.5 4],...
        'paperpositionmode','auto');
    
        tiledlayout('flow')
        nexttile
        boxplot(g1_cc,labels)
        ylabel('Nodal Clustering Coefficient')
        title({' ',g1t1_struct.grp})
        
        nexttile
        boxplot(g2_cc,labels)
        title({' ',g2t1_struct.grp})
        
        sgtitle({[label_ci ' ref-partition w/in module nodal clustering coefficient'],['network thresh p<' num2str(pcmp) ' module ' num2str(ciu(cmi))]},'Interpreter','none')
        
        clear g1_cc g2_cc
        
        filename = fullfile(outdir, [label_ci '_refPart_within_clustCoef_' g1t1_struct.grp '_vs_' g2t1_struct.grp '_thr' num2str(pcmp) '_m' num2str(ciu(cmi)) '.pdf']);
        print(gcf,'-dpdf',filename)
        clear filename
        
    end
    clear pcmp
    
end
close all
end

