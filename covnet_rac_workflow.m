%% COPY THIS SCRIPT ONE LEVEL UP FROM THE COVNET DIRECTORY AND NAME IT SO IT MAKES SENSE FOR YOUR PROJECT.

%% YOU SHOULD HAVE A PROJECT DIRECRORY WITH THE COVNET CODE FOLDER AND THE WORKFLOW SCRIPT WHEN YOU START.

%% THEN FOLLOW THE STEP BELOW TO CREATE AND SAVE THE DATA CELLS AS .MAT FILES IN THE PROJECT DIRECTORY.

%% THE PROJECT DIRECTORY SHOULD BE YOUR WORKING DIRECTORY WHEN YOU RUN THE SCRIPT, SO THAT ALL OUTPUTS ARE
%% WRITTEN EITHER TO THE PROJECT DIRECTORY OR TO SUBDIRECTORIES CREATED BY THIS CODE.


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
addpath(genpath('CovNet'))

    % Load region labels for your data.
%load('roi_labels.mat')   
    
    % This will load colormaps contained in the package
load colormaps.mat

    % Thresholding is an option for covariance networks. P-value thresholds are
    % set here and permutation testing is done to retain only edges with p-values
    % below threshold. A single, multiple, or no p-values can be provided.
pval = [0.05];

%% Run tier 1 analysis
%   - generates covariance matrices (thresholded and unthresholded)
%   - runs multiresolution modularity
%   - computes adjusted mutual information among modularity partitions
%   - pairwise permutation testing (work in progress)
%   - computes global network metrics

% load and run data
load('PPE_RAC.mat')
load('ppe_covariates_n39.mat','PPE_COVS3')
covariance_analysis_tier1(PPE_RAC,roi,bluered_cmap,pval,PPE_COVS3)

% running a 3 group variant with nonsmoking CON and smoking and nonsmoking
% NTS
load('PPE_RAC_3smkGRP_age_sex.mat')
covariance_analysis_tier1(PPE_3smkGRP,roi,bluered_cmap,pval,PPE_3smkGRP_COVS)

% load and run data
load('PPE_roiSize.mat')
load('ppe_covariates_n39.mat','PPE_COVS3')
covariance_analysis_tier1(PPE_roiSize,roi_labels,bluered_cmap,pval,PPE_COVS3)

%% Extract summary metrics from modules for a given reference partition
    % Set your reference modularity partition from a tier 1 output.
    
    % load tier1 data for a gorup    
ctrl_tier1 = load('covariance_out_tier1_PPE_age_smk_sex.mat');
ctrl_3grp_tier1 = load('covariance_out_tier1_PPE_3grp_age_sex.mat');
    
    % Set Controls as reference partition 
ci = ctrl_tier1.mrccPartition{2,2};
ci_3grp=ctrl_3grp_tier1.mrccPartition{2,2};

    % Set the groups from which you want to extract module specific summaries
    % example only the males
cmask = [1;1];
cmask_3grp = [1;1;1];
    % The variable will look like this: 
    % (notice the size is the same as the size of data entries in the cell structure)
    % 1 0
    % 1 0
    % 1 0
    % it should be the same size as number of data cells (in our example
    % its 3 ages on rows and 2 sexes on columns).

    % Set output file root for a text output that is written out:
outroot1 = 'ref_ctrl';
outroot2 = 'ref_3grp_ctrl';

        % Run the function
    % The function extracts mean regional values for specified reference modules.
    % Mean and StDev for each module are saved to the text file.
fcn_summary_modules(ci,ctrl_tier1,cmask,outroot1)

fcn_summary_modules(ci_3grp,ctrl_3grp_tier1,cmask_3grp,outroot2)
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
covariance_analysis_tier2(ctrl_tier1,xfad_tier1,cmask)

    












