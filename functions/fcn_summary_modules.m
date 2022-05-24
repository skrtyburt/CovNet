function fcn_summary_modules(ci,tier1_out_group1,cmask,outprefix)

narginchk(2,4)

% check comm structure
 if ~isstruct(tier1_out_group1)
    fprintf(2,'tier1_out is not a structure. Exiting...')
    return
 end
% check node size vs size of community structure vector
if tier1_out_group1.N ~= length(ci)
    fprintf(2,'number of nodes does not match between inputs. Exiting...')
    return
end
% make sure communities are linearly indexed
if sum(ci==1)==0 || length(unique(diff(unique(ci)))) > 1 
    % if index doesnt start with one or if there are gaps in index
    ci = fcn_relabel_partitions(ci);
end
% optional cmask should be the same size as number of DATA cells in the
% original input cell structure
if exist('cmask','var')
    [r,c]=size(cmask);
    if r ~= tier1_out_group1.Nr-1 || c ~= tier1_out_group1.Nc-1
        fprintf(2,'size of cmask does not match size data cells')
        return
    end
else
    cmask=ones(tier1_out_group1.Nr-1,tier1_out_group1.Nc-1);
end
% find cells from which to compute summary values
cmask = vertcat(zeros(1,tier1_out_group1.Nc), horzcat(zeros(tier1_out_group1.Nr-1,1),cmask));
[rw,cl]=find(cmask);
% find unique community labels
ciu = unique(ci);

% find largest group size
slength = max(max(tier1_out_group1.S));

% output
gmeans = table;
% for every indexed group
for i=1:length(rw)
    d = tier1_out_group1.cellData{rw(i),cl(i)};
    % for every community label
    for cj=1:length(ciu)
        % extract the average regional value in community
        cmean = mean(d(ci==ciu(cj),:));
        cstd = std(d(ci==ciu(cj),:));
        nm=[tier1_out_group1.cellData{rw(i),1} '_' tier1_out_group1.cellData{1,cl(i)} '_module' num2str(ciu(cj))];
        if length(cmean)~=slength
            cmean(end+1:slength)=nan;
            cstd(end+1:slength)=nan;
        end
        nm1=['mn_' nm];
        gmeans.(nm1)=cmean';
        nm2=['std_' nm];
        gmeans.(nm2)=cstd';
        clear cmean nm nm1 nm2
    end
    clear d Md
end

if exist('outprefix','var')
    writetable(gmeans,[outprefix '_comsummary.txt'],'Delimiter','\t')
else
    writetable(gmeans,'refMod_comsummary.txt','Delimiter','\t')
end












