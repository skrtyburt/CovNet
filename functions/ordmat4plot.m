


function [boxData, xaxisLabel] = ordmat4plot(m1,m2,cNames,gNames)

r = length(cNames)*2;

for i =1:length(cNames)
    ii = 1:2:r;
    boxData(ii(i),:) = m1(:,i)';
    xaxisLabel{ii(i),1} = [gNames{1} '-' cNames{i}];
    jj = 2:2:r;
    boxData(jj(i),:) = m2(:,i)';
    xaxisLabel{jj(i),1} = [gNames{2} '-' cNames{i}];
end

