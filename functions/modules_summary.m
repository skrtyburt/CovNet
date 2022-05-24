function modules_summary(data, partition)

narginchk(2,2)
[N,S,G] = size(data);

if N~=length(partition)
    fprintf(2,'Number of nodes in parition does not match data.\n')
    return
end
if ~ismatrix(partition) || min(size(partition))>1
    fprintf(2,'Partition must be a [N,1] size vector.\n')
    return
end

%%% make boxplots of within module averages of network properties and save
%%% data

%%% compute participation coeffcient and compare between networks and save

%%% plot mean regional suvR for modules between groups and save in a stats
%%% friendly format