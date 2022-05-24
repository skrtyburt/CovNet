function fig_scatter_by_group(data_vector,grouping_vector,labels)
% FIG_SCATTER_BY_GROUP      plot individual points bar graph style
%                           with a line at mean
%
%   fig_scatter_by_group(data,groups,labels)
%
%   Inputs:
%       data_vector - column of values, where each subject is represented
%                       once.
%
%       grouping_vector - column of integer values designating group
%                           membership.
%
%       labels - a cell with text lables for the groups ordered
%                   corresponding to increasing values in grouping_vector.
%
%   Evgeny Jenya Chumin, Indiana University 2020
%
%   Modification History:
%       2020: Original
%%
uG = unique(grouping_vector);                       % groups
for g=1:length(uG)
    numS = sum(grouping_vector==uG(g));             % subjects in group
    
    x = g-.25:((.25*2)/numS):g+.25;       % distribute points on the x-axis
    x = x(1:end-1)';
    
    scatter(x,data_vector(grouping_vector==uG(g)))  % plot
    hold on
    
    %draw a line at group mean
    mean_line = zeros(numS,1)+nanmean(data_vector(grouping_vector==uG(g)));
    plot(x,mean_line,'k','LineWidth',2)
end

xticks(1:1:length(uG))      % label the groups with input
if exist('labels','var')
    xticklabels(labels)     % use input cell of labels
else
    xticklabels(uG);        % use unique integers of grouping_vector
end
xlim([0 length(uG)+1])