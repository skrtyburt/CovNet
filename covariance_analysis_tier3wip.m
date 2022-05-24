function covariance_analysis_tier3(group1_tier1out,group2_tier1out)

narginchk(2,2)

g1 = group1_tier1out;
g2 = group2_tier1out;

if g1.N ~= g2.N
    fprintf(2,'Unequal number of nodes in group networks. Exiting ...\n')
    return
end
N=g1.N;
if g1.Nr ~= g2.Nr
    fprintf(2,'Unequal number of row subgroups in cell. Exiting ...\n')
    return
elseif g1.Nc ~= g2.Nc
    fprintf(2,'Unequal number of column subgroups in cell. Exiting ...\n')
    return
end
Nr=g1.Nr; 
Nc=g1.Nc;
ccomp = unique(strcmp(g1.cNames(2:end),g2.cNames(2:end)));
rcomp = unique(strcmp(g1.rNames(2:end),g2.rNames(2:end)));
if length(ccomp)>1 || length(rcomp)>1
    fprintf(2,'Labels in group cells do not match. Exiting ...\n')
    return
end
rNames = g1.rNames(2:end); 
cNames = g1.cNames(2:end);
clear ccomp rcomp

gNames{1}=g1.grp;
gNames{2}=g2.grp;


c1 = g1.covMat;
c2 = g2.covMat;


% node strenth (unthresholded)

g1.Spos = c1(1:end,1);  g1.Spos(1,2:Nc) = cNames;  
g1.Sneg = c1(1:end,1);  g1.Sneg(1,2:Nc) = cNames;

g2.Spos = c2(1:end,1);  g2.Spos(1,2:Nc) = cNames;  
g2.Sneg = c2(1:end,1);  g2.Sneg(1,2:Nc) = cNames;

for i = 2:size(c1,1)
    for j = 2:size(c1,2)
        [g1.Spos{i,j},g1.Sneg{i,j},g1_totalSpos(i-1,j-1),g1_totalSneg(i-1,j-1)]=strengths_und_sign(c1{i,j});
        [g2.Spos{i,j},g2.Sneg{i,j},g2_totalSpos(i-1,j-1),g2_totalSneg(i-1,j-1)]=strengths_und_sign(c2{i,j});
    end
end

% bar plot total strength
%positive
[bd1,xL] = ordmat4plot(g1_totalSpos,g2_totalSpos,cNames,gNames);
[bd2,~] = ordmat4plot(g1_totalSneg,g2_totalSneg,cNames,gNames);

figure
subplot(1,2,1)
bar(bd1); xticklabels(xL); ylabel('Sum Positive Nodal Strength')
title('Positive Strength')
legend(rNames,'Location','northeastoutside')
subplot(1,2,2)
bar(bd2); xticklabels(xL); ylabel('Sum Negative Nodal Strength')
title('Negative Strength')
legend(rNames,'Location','northeastoutside')



























