%% Edgewise Fisher testing for significant difference between correlations
FisherP = cell.empty;
ci = 1; % counting indices
for r=2:Nr % every row (g1)
for c=2:Nc % every column (g1)
    ci=ci+1;
    FisherP{ci,1}=[grp '_' cellData{r,1} '_' cellData{1,c}];
    cii = 1;
    for rr=2:Nr % every row (g2)
    for cc=2:Nc % every column (g2)
        cii=cii+1;
        FisherP{1,cii}=[grp '_' cellData{rr,1} '_' cellData{1,cc}];
        FisherP{ci,cii}=double.empty;
        r1 = covMat{r,c};           % set covariance matrix for group1
        n1 = S(r-1,c-1);            % set number of animals in group1
        r2 = covMat{rr,cc};         % covariance group2
        n2 = S(rr-1,cc-1);          % number of animals group2
        for ii=1:length(r1) % loop over rows
            for jj=1:length(r1) % loop over columns
                % test and save p-values
                [FisherP{ci,cii}(ii,jj)] = fisher_r2z_compare(r1(ii,jj),n1,r2(ii,jj),n2);
            end
        end
        clear r1 n1 r2 n2
    end
    end
end
end