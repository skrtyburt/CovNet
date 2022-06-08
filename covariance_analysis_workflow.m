
%% Step1, Prepare data:
% example:
% control animals at 3 ages and two sexes.
% Currently organized as animal by region.

% MAKE SURE THAT DATA THAT GOES INTO THE CELL ARE 
% REGION (ROWS) BY ANIMAL (COLUMN).
%%
% Example of how cell structure was buld:
%{

 1. Load data that was imported from excel
    each group was pasted into a separate variable in matlab and save to
    the file that is being read in

load('CTRL.mat') 

 2. Organize data into a cell:

    % Initialize cell
CTRL = cell.empty;

    % Group title into the first row/column of the cell
CTRL{1,1}='Control';

    % Assign factor labels to the first column (here age)
CTRL{2,1}='4m';
CTRL{3,1}='6m';
CTRL{4,1}='12m';

    % Assign factor labels to the first column (here sex)
CTRL{1,2}='M';
CTRL{1,3}='F';

    % Assign the data into the corresponding {row,column} indices
        % Note: HERE I AM TRANSPOSING THE DATA (using ' ) SO THAT ROWS ARE REGIONS
CTRL{2,2}=CTRL_4M';  
CTRL{2,3}=CTRL_4F';

CTRL{3,2}=CTRL_6M';  
CTRL{3,3}=CTRL_6F';

CTRL{4,2}=CTRL_12M'; 
CTRL{4,3}=CTRL_12F';

 3. Save the cell structure
   
    % save to a new file
save('Control_data.mat','CTRL')

    % append an existing file (e.g. store the cell structure in the same
    file as source data.
save('CTRL.mat','CTRL','-append')

%}

%% Prepare the necessary inputs

    % Add code package to matlab path. (THIS PATH WILL VARY FOR EACH USER)
addpath(genpath('../CovNet'))

    % Load region labels for your data.
load('roi_labels.mat')   
    
    % This will load colormaps contained in the package
load colormaps.mat

    % Thresholding is an option for covariance networks. P-value thresholds are
    % set here and permutation testing is done to retain only edges with p-values
    % below threshold. A single, multiple, or no p-values can be provided.
pval = [0.05 0.01];

%% Run tier 1 analysis
%   - generates covariance matrices (thresholded and unthresholded)
%   - runs multiresolution modularity
%   - computes adjusted mutual information among modularity partitions
%   - pairwise permutation testing (work in progress)
%   - computes global network metrics

% load and run control data
load('CTRL.mat')
covariance_analysis_tier1(CTRL,roi_labels,bluered_cmap,pval)

% load and run mutant data
load('XFAD.mat')
covariance_analysis_tier1(XFAD,roi_labels,bluered_cmap,pval)

%% Extract summary metrics from modules for a given reference partition
    % Set your reference modularity partition from a tier 1 output.
    
    % load tier1 data for a gorup    
ctrl_tier1 = load('covariance_out_tier1_Control.mat');
    
    % Set as reference the 4month Male Control Consensus partition
ci = ctrl_tier1.mrccPartition{2,2};

    % Set the groups from which you want to extract module specific summaries
    % example only the males
cmask = [1,0; 1,0; 1,0];
    % The variable will look like this: 
    % (notice the size is the same as the size of data entries in the cell structure)
    % 1 0
    % 1 0
    % 1 0
    % it should be the same size as number of data cells (in our example
    % its 3 ages on rows and 2 sexes on columns).

    % Set output file root for a text output that is written out:
outroot1 = 'ctrl_4M_M_ref_ctrl';

        % Run the function
    % The function extracts mean regional values for specified reference modules.
    % Mean and StDev for each module are saved to the text file.
fcn_summary_modules(ci,ctrl_tier1,cmask,outroot1)
%-------------------------------------------------------------------------%
    % Repeat for XFAD using the 4month Male XFAD as reference partition and the same cmask:
xfad_tier1 = load('covariance_out_tier1_xfad.mat');
    % if you want to change reference to xfad 4m M
ci2 = xfad_tier1.mrccPartition{2,2};
    % label for the text file
outroot2 = 'xfad_4M_M_ref_xfad';
    % run the script
fcn_summary_modules(ci2,xfad_tier1,cmask,outroot2)

%% tier 2 compare network metrics

     % Compare Males at each age (Set by the cmask input) between control
     % and xfad at p-value thresholds set in tier1.
     
     % Note: Network metric comparisons are only done on thresholded networks. 
covariance_analysis_tier2a(ctrl_tier1,xfad_tier1,cmask)

    % Compare within module network metrics (given a reference partition)













