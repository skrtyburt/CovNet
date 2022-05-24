function [P, Q] = louvain_best(mat,gamma_vector,iterations)
% LOUVAIN_BEST      select highest Q partition from iterations
%
% Requirements:     bct_community_louvain.m
%
% Notes: 
%   Defaults to negative assymetry model null model for modularity.
%
%   [P,Q] = louvain_best(mat,gamma_vector,iterations)
%
%   Inputs:
%       mat -   NxN adjacency matrix 
%               (tested and ran on functional matrices)
%       
%       gamma_vector -  resolution parameter values to be ran
%                       all gamma must be y>0
%
%       iterations -    number of repeats of louvain algorithm at each
%                       gamma.
%
%   Outputs:
%       P - community affilication vectors.
%           a matrix size (N,Y), where N is the number of nodes in mat
%           anf Y is the number of gamma values ran.
%
%       Q - modularity statistic.
%           a matrix size (N,Y), where N is the number of nodes in mat
%           anf Y is the number of gamma values ran.
%
%   Example Inputs:
%       mat is the input connectivity matrix.
%       gamma_vector = 0.1:.1:4;    % range of gamma resolution values
%                                    % starting from .1 to 4 in .1 increments
%       iterations = 1000;          % number of repeats at each resolution
%
%   Evgeny Jenya Chumin, Indiana University 2019
%
%   Modification History:
%       2019: Original
%%
N = size(mat,1);
%
idx=sort(find(nansum(mat)==0),'descend');           % find any NaN nodes
%   
for ii=1:length(idx)                                % remove NaN rows from matrix
    mat=vertcat(mat(1:idx(ii)-1,:),mat(idx(ii)+1:end,:));
end
for ii=1:length(idx)                                % remove NaN columns form matrix
    mat=horzcat(mat(:,1:idx(ii)-1),mat(:,idx(ii)+1:end));
end
%
numY = length(gamma_vector);
Nv = size(mat,1);
%
mat(eye(Nv,Nv)==1)=0;                               % zero diagonal
idx=sort(idx);                                      % reverse sort order for adding NaNs later
%%

%
Mr=nan(Nv,iterations);
Qr=nan(iterations,1);
Q=nan(1,numY);
P=nan(N,numY);
%
for y=1:numY
    res = gamma_vector(y);
    fprintf('gamma: %d/%d\n',y,numY)
    parfor r=1:iterations
        [Mr(:,r),Qr(r,1)] = community_louvain(mat,res,[],'negative_asym');
    end
    Qmax=max(Qr);
    Q(1,y)=Qmax;
    Qidx=find(Qr==Qmax);
    if length(Qidx)>1                               % if max Q is for more than 1 partition
        for p=1:length(Qidx)                        % isolate maxQ partitions
            Mfi(:,p)=Mr(:,Qidx(p));
        end
        [uMfi,~,ic]=unique(Mfi','rows','stable');   % find unique partitions
        clear Mfi
        counts = accumarray(ic,1);                  % count the unique partitions
        bpart=find(counts==max(counts));
        if length(bpart) > 1                        % if more than 1 maxQ partition occurs with equal frequency
            bpart=bpart(1);                         % take the first
        end
        Mb=uMfi(bpart,:)';                           % take the most frequent maxQ
        clear uMfi
    elseif length(Qidx)==1                          % if only one maxQ partition, take that one
        Mb=Mr(:,Qidx);
    end
    for ii=1:length(idx)                            % add back any NaN nodes
        Mb=vertcat(Mb(1:idx(ii)-1,1),nan(1,1),Mb(idx(ii):end,1));
    end
    P(:,y)=Mb;
    clear Mb Mr uMfi Mfi Qidx Qmax
end 
 delete(gcp('nocreate'))