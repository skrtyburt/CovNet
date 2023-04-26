%% YOU SHOULD HAVE A PROJECT DIRECRORY WITH THE COVNET CODE FOLDER AND THE WORKFLOW SCRIPT WHEN YOU START.
%% COPY THIS SCRIPT ONE LEVEL UP FROM THE MAIN COVNET DIRECTORY AND NAME IT SO IT MAKES SENSE FOR YOUR PROJECT.

%% IF YOU DATA IS NOT IN MATLAB OR NEEDS TO BE FORMATTER:
% FOLLOW THE STEP BELOW TO CREATE AND SAVE THE DATA CELLS AS .MAT FILES IN THE PROJECT DIRECTORY.
%% THE PROJECT DIRECTORY SHOULD BE YOUR WORKING DIRECTORY WHEN YOU RUN THE SCRIPT, SO THAT ALL OUTPUTS ARE
%%     WRITTEN EITHER TO THE PROJECT DIRECTORY OR TO SUBDIRECTORIES CREATED BY THIS CODE.
%---------------------------------------------------------------------------------------------------------%
%% Step 1, Prepare data:
    % Example: control animals at 3 ages and two sexes.
    % Currently organized as animal (rows) by region (columns).

%% MAKE SURE THAT DATA THAT GOES INTO THE CELL ARE REGION (ROWS) BY ANIMAL (COLUMN)!

%% Example of how cell structure was buld:
    %{
    % 1. Load data that was imported from excel. Each group was pasted into a separate variable in matlab
    %    and save to the file that is being read in
         load('CTRL.mat') 
    % 2. Organize data into a cell:
    % Initialize cell
         CTRL = cell.empty;
    % Group or study title into the first row/column of the cell
        CTRL{1,1}='Control';
    % Assign labels in the first column (here age)
        CTRL{2,1}='4m';
        CTRL{3,1}='6m';
        CTRL{4,1}='12m';
    % Assign labels in the first row (here sex)
        CTRL{1,2}='M';
        CTRL{1,3}='F';
    % Assign the data into the corresponding {row,column} indices
    % NOTE: Here I am transposing the data (using ' ) so that rows are
    %       regions and columns are subjects/animals.
        CTRL{2,2}=CTRL_4M';  
        CTRL{2,3}=CTRL_4F';
        CTRL{3,2}=CTRL_6M';  
        CTRL{3,3}=CTRL_6F';
        CTRL{4,2}=CTRL_12M'; 
        CTRL{4,3}=CTRL_12F';
    % 3. Save the cell structure
    % save to a new file 
        save('Control_data.mat','CTRL')  % OR
    % append an existing file (e.g. store the cell structure in the same
    % file as source data.
        save('CTRL.mat','CTRL','-append')
    %}

%% Prepare the necessary inputs
% Add code package to matlab path. (THIS PATH WILL VARY FOR EACH USER)
    addpath(genpath('C:\Users\echumin\Documents\GitHub\CovNet'))

% Load region labels for your data.
  %  load('roi_labels.mat')  % a cell of strings
  % Load V2 which has groupings of nodes into systems/networks
        load('roiLabels_wGroupings.mat')
    
% This will load colormaps contained in the package
    load colormaps.mat

% Thresholding is an option for covariance networks. P-value thresholds are
% set here and permutation testing is done to retain only edges with p-values
% below threshold. A single, multiple, or no p-values can be provided.
% pval =[];
  pval = [0.05 0.01];

%% Run tier 1 analysis
% Tier 1 is ran on a single cellData. If you have multiple cellData
% elements that you want to compare (e.g. here there is a CTRL and 5XFAD
% genotype groups, each with 2 sexes at 3 ages), they are ran separetly and
% can be compared in Tier 2.

    %   - generates covariance matrices (thresholded and unthresholded)
    %   - runs multiresolution modularity - MRCC - Jeub, Sporns, Fortunat0, 2018
    %   - computes adjusted mutual information among modularity partitions
    %   - computes global network metrics
    %   - group difference testing

% load and run control data
load('CTRL.mat')
covariance_analysis_tier1(CTRL,roi_labels,bluered_cmap,pval)

% load and run mutant data
load('XFAD.mat')
covariance_analysis_tier1(XFAD,roi_labels,bluered_cmap,pval)

%% Extract summary metrics from modules for a given reference partition
    
% Set your reference modularity partition from a tier 1 output.
    % CTRL %
        ctrl_tier1 = load('covariance_out_tier1_PNG_r600_Control.mat');
    % XFAD %
        xfad_tier1 = load('covariance_out_tier1_PNG_r600_xfad.mat');
        
%% perform statistics on module average data
% CTRL
% 4month MALES (REF = CTRL 4mo MALE)
    % Set the groups from which you want to extract module specific summaries 
        % only the males
            cmask_all = [1,1; 1,1; 1,1];
            % this variable will look like this:
            % 1 1
            % 1 1
            % 1 1
            % it should be the same size as number of data cells (in our example
            % its 3 ages on rows and sex as two columns).
    % set output file root for any outputs that are written out
            % first example below
        outroot1 = '1_CTRL_ref_ctrl_4mo_M';
    % set the reference partition
        ci_ctrl_4mo_M = ctrl_tier1.mrccPartition{2,2};
    % run the script
    % function for extracting mean SUVR for specified input modules
    % means and stdev for each module are saved to a text file
        fcn_summary_modules(ci_ctrl_4mo_M,ctrl_tier1,cmask_all,outroot1)

        outroot2 = '2_XFAD_ref_ctrl_4mo_M';
        ci_xfad_4mo_M = xfad_tier1.mrccPartition{2,2};
        fcn_summary_modules(ci_xfad_4mo_M,xfad_tier1,cmask_all,outroot2)


        outroot3 = '3_joinedCTRL_XFAD_ref_ctrl_4mo_M';
        fcn_summary_modules(ci_ctrl_4mo_M,ctrl_tier1,cmask_all,outroot3,xfad_tier1)

% Exploratory: get system average SUVr with different system groupings
        outroot4 = '4_macnicol';
        fcn_summary_modules(macnicolGRPs',ctrl_tier1,cmask_all,outroot4,xfad_tier1,1)
        outroot5 = '5_func10';
        fcn_summary_modules(functionalGRPs',ctrl_tier1,cmask_all,outroot5,xfad_tier1,1)
        outroot6 = '6_domain7';
        fcn_summary_modules(domainGRPs',ctrl_tier1,cmask_all,outroot6,xfad_tier1,1)



        

        
% (2) all CTRL FEMALES (REF = CTRL 4mo FEMALE)
    % Set the groups from which you want to extract module specific summaries 
        % only the females
            cmask_F = [0,1; 0,1; 0,1];
    % set output file root for any outputs that are written out
        outroot2 = '2_CTRL_FEMALES_ref_ctrl_4mo_F';
    % set the reference partition
        ci_ctrl_4mo_M = ctrl_tier1.mrccPartition{2,3};
    % run the script
        fcn_summary_modules(ci_ctrl_4mo_M,ctrl_tier1,cmask_F,outroot2)

% XFAD
% (3) all XFAD MALES (REF = XFAD 4mo MALE)
   % Set the groups from which you want to extract module specific summaries 
        % only the males
            cmask_M = [1,0; 1,0; 1,0];
    % set output file root for any outputs that are written out
        outroot3 = '3_XFAD_MALES_ref_xfad_4mo_M';
    % set the reference partition
        ci_xfad_4mo_M = xfad_tier1.mrccPartition{2,2};
    % run the script
        fcn_summary_modules(ci_xfad_4mo_M,xfad_tier1,cmask_M,outroot3)       % fixed typo: uppercae M

% (4) all XFAD FEMALES (REF = XFAD 4mo FEMALE)
    % Set the groups from which you want to extract module specific summaries 
        % only the females
            cmask_F = [0,1; 0,1; 0,1];
    % set output file root for any outputs that are written out
        outroot4 = '4_XFAD_FEMALES_ref_xfad_4mo_F';
    % set the reference partition
        ci_xfad_4mo_F = xfad_tier1.mrccPartition{2,3};
    % run the script
        fcn_summary_modules(ci_xfad_4mo_F,xfad_tier1,cmask_F,outroot4)

        
%% Genotype model(s)        n=6 (5-10)

% MALE
% (5) XFAD 4mo MALE (REF = CTRL 4mo MALE)
    % Set the groups from which you want to extract module specific summaries
        % only 4mo MALES
            cmask_4mo_M = [1,0; 0,0; 0,0];
    % set output file root for any outputs that are written out
        outroot5a = '5a_XFAD_4MO_MALE_ref_ctrl_4mo_M';
    % set the reference partition
        ci_ctrl_4mo_M = ctrl_tier1.mrccPartition{2,2};
    % run the script
        fcn_summary_modules(ci_ctrl_4mo_M,xfad_tier1,cmask_4mo_M,outroot5a)   
    % set output file root for any outputs that are written out
        outroot5b = '5b_CTRL_4MO_MALE_ref_ctrl_4mo_M';
    % set the reference partition
        ci_ctrl_4mo_M = ctrl_tier1.mrccPartition{2,2};
    % run the script
        fcn_summary_modules(ci_ctrl_4mo_M,ctrl_tier1,cmask_4mo_M,outroot5b) 
%-----------------------------------------%
ci_xfad_4mo_M = xfad_tier1.mrccPartition{2,2};
cmask_4mo_M = [1,0; 0,0; 0,0];
outroot5a2 = 'xfad_ref_4genotype/5a2_XFAD_4MO_MALE_ref_xfad_4mo_M';
fcn_summary_modules(ci_xfad_4mo_M,xfad_tier1,cmask_4mo_M,outroot5a2)   

outroot5b2 = 'xfad_ref_4genotype/5b2_CTRL_4MO_MALE_ref_xfad_4mo_M';
fcn_summary_modules(ci_xfad_4mo_M,ctrl_tier1,cmask_4mo_M,outroot5b2)
%-----------------------------------------%

% (6) XFAD 6mo MALE (REF = CTRL 6mo MALE)
    % Set the groups from which you want to extract module specific summaries
        % only 6mo MALES
            cmask_6mo_M = [0,0; 1,0; 0,0];
    % set output file root for any outputs that are written out
        outroot6a = '6a_XFAD_6MO_MALE_ref_ctrl_6mo_M';
    % set the reference partition
        ci_ctrl_6mo_M = ctrl_tier1.mrccPartition{3,2};
    % run the script
        fcn_summary_modules(ci_ctrl_6mo_M,xfad_tier1,cmask_6mo_M,outroot6a) 
    % set output file root for any outputs that are written out
        outroot6b = '6b_CTRL_6MO_MALE_ref_ctrl_6mo_M';
    % set the reference partition
        ci_ctrl_6mo_M = ctrl_tier1.mrccPartition{3,2};
    % run the script
        fcn_summary_modules(ci_ctrl_6mo_M,ctrl_tier1,cmask_6mo_M,outroot6b)
        
%-----------------------------------------%
ci_xfad_6mo_M = xfad_tier1.mrccPartition{3,2};
cmask_6mo_M = [0,0; 1,0; 0,0];
outroot6a2 = 'xfad_ref_4genotype/6a2_XFAD_6MO_MALE_ref_xfad_6mo_M';
fcn_summary_modules(ci_xfad_6mo_M,xfad_tier1,cmask_6mo_M,outroot6a2)   

outroot6b2 = 'xfad_ref_4genotype/6b2_CTRL_6MO_MALE_ref_xfad_6mo_M';
fcn_summary_modules(ci_xfad_6mo_M,ctrl_tier1,cmask_6mo_M,outroot6b2)
%-----------------------------------------%
        
% (7) XFAD 12mo MALE (REF = CTRL 12mo MALE)
    % Set the groups from which you want to extract module specific summaries
        % only 12mo MALES
            cmask_12mo_M = [0,0; 0,0; 1,0];
    % set output file root for any outputs that are written out
        outroot7a = '7a_XFAD_12MO_MALE_ref_ctrl_12mo_M';
    % set the reference partition
        ci_ctrl_12mo_M = ctrl_tier1.mrccPartition{4,2};
    % run the script
        fcn_summary_modules(ci_ctrl_12mo_M,xfad_tier1,cmask_12mo_M,outroot7a)  
    % set output file root for any outputs that are written out
        outroot7b = '7b_CTRL_12MO_MALE_ref_ctrl_12mo_M';
    % set the reference partition
        ci_ctrl_12mo_M = ctrl_tier1.mrccPartition{4,2};
    % run the script
        fcn_summary_modules(ci_ctrl_12mo_M,ctrl_tier1,cmask_12mo_M,outroot7b)

%-----------------------------------------%
ci_xfad_12mo_M = xfad_tier1.mrccPartition{4,2};
cmask_12mo_M = [0,0; 0,0; 1,0];
outroot7a2 = 'xfad_ref_4genotype/7a2_XFAD_12MO_MALE_ref_xfad_12mo_M';
fcn_summary_modules(ci_xfad_12mo_M,xfad_tier1,cmask_12mo_M,outroot7a2)   

outroot7b2 = 'xfad_ref_4genotype/7b2_CTRL_12MO_MALE_ref_xfad_12mo_M';
fcn_summary_modules(ci_xfad_12mo_M,ctrl_tier1,cmask_12mo_M,outroot7b2)
%-----------------------------------------%

% FEMALE
% (8) XFAD 4mo FEMALE (REF = CTRL 4mo FEMALE)
    % Set the groups from which you want to extract module specific summaries
        % only 4mo FEMALES
            cmask_4mo_F = [0,1; 0,0; 0,0];
    % set output file root for any outputs that are written out
        outroot8a = '8a_XFAD_4MO_FEMALE_ref_ctrl_4mo_F';
    % set the reference partition
        ci_ctrl_4mo_M = ctrl_tier1.mrccPartition{2,3};
    % run the script
        fcn_summary_modules(ci_ctrl_4mo_M,xfad_tier1,cmask_4mo_F,outroot8a)  
    % set output file root for any outputs that are written out
        outroot8b = '8b_CTRL_4MO_FEMALE_ref_ctrl_4mo_F';
    % set the reference partition
        ci_ctrl_4mo_M = ctrl_tier1.mrccPartition{2,3};
    % run the script
        fcn_summary_modules(ci_ctrl_4mo_M,ctrl_tier1,cmask_4mo_F,outroot8b)

%-----------------------------------------%
ci_xfad_4mo_F = xfad_tier1.mrccPartition{2,3};
cmask_4mo_F = [0,1; 0,0; 0,0];
outroot8a2 = 'xfad_ref_4genotype/8a2_XFAD_4MO_FEMALE_ref_xfad_4mo_F';
fcn_summary_modules(ci_xfad_4mo_F,xfad_tier1,cmask_4mo_F,outroot8a2)   

outroot8b2 = 'xfad_ref_4genotype/8b2_CTRL_4MO_FEMALE_ref_xfad_4mo_F';
fcn_summary_modules(ci_xfad_4mo_F,ctrl_tier1,cmask_4mo_F,outroot8b2)
%-----------------------------------------%

% (9) XFAD 6mo FEMALE (REF = CTRL 6mo FEMALE)
    % Set the groups from which you want to extract module specific summaries
        % only 6mo FEMALES
            cmask_6mo_F = [0,0; 0,1; 0,0];
    % set output file root for any outputs that are written out
        outroot9a = '9a_XFAD_6MO_FEMALE_ref_ctrl_6mo_F';
    % set the reference partition
        ci_ctrl_6mo_F = ctrl_tier1.mrccPartition{3,3};
    % run the script
        fcn_summary_modules(ci_ctrl_6mo_F,xfad_tier1,cmask_6mo_F,outroot9a)
    % set output file root for any outputs that are written out
        outroot9b = '9b_CTRL_6MO_FEMALE_ref_ctrl_6mo_F';
    % set the reference partition
        ci_ctrl_6mo_F = ctrl_tier1.mrccPartition{3,3};
    % run the script
        fcn_summary_modules(ci_ctrl_6mo_F,ctrl_tier1,cmask_6mo_F,outroot9b)

%-----------------------------------------%
ci_xfad_6mo_F = xfad_tier1.mrccPartition{3,3};
cmask_6mo_F = [0,0; 0,1; 0,0];
outroot9a2 = 'xfad_ref_4genotype/9a2_XFAD_6MO_FEMALE_ref_xfad_6mo_F';
fcn_summary_modules(ci_xfad_6mo_F,xfad_tier1,cmask_6mo_F,outroot9a2)   

outroot9b2 = 'xfad_ref_4genotype/9b2_CTRL_6MO_FEMALE_ref_xfad_6mo_F';
fcn_summary_modules(ci_xfad_6mo_F,ctrl_tier1,cmask_6mo_F,outroot9b2)
%-----------------------------------------%

% (10) XFAD 12mo FEMALE (REF = CTRL 12mo FEMALE)
    % Set the groups from which you want to extract module specific summaries
        % only 12mo FEMALES
            cmask_12mo_F = [0,0; 0,0; 0,1];
    % set output file root for any outputs that are written out
        outroot10a = '10a_XFAD_12MO_FEMALE_ref_ctrl_12mo_F';
    % set the reference partition
        ci_ctrl_12mo_F = ctrl_tier1.mrccPartition{4,3};
    % run the script
        fcn_summary_modules(ci_ctrl_12mo_F,xfad_tier1,cmask_12mo_F,outroot10a)
    % set output file root for any outputs that are written out
        outroot10b = '10b_CTRL_12MO_FEMALE_ref_ctrl_12mo_F';
    % set the reference partition
        ci_ctrl_12mo_F = ctrl_tier1.mrccPartition{4,3};
    % run the script
        fcn_summary_modules(ci_ctrl_12mo_F,ctrl_tier1,cmask_12mo_F,outroot10b)

%-----------------------------------------%
ci_xfad_12mo_F = xfad_tier1.mrccPartition{4,3};
cmask_12mo_F = [0,0; 0,0; 0,1];
outroot10a2 = 'xfad_ref_4genotype/10a2_XFAD_12MO_FEMALE_ref_xfad_12mo_F';
fcn_summary_modules(ci_xfad_12mo_F,xfad_tier1,cmask_12mo_F,outroot10a2)   

outroot10b2 = 'xfad_ref_4genotype/10b2_CTRL_12MO_FEMALE_ref_xfad_12mo_F';
fcn_summary_modules(ci_xfad_12mo_F,ctrl_tier1,cmask_12mo_F,outroot10b2)
%-----------------------------------------%

        
%% Sex model(s)     n=6 (11-16)

% CTRL
% (11) CTRL 4mo FEMALE (REF = CTRL 4mo MALE)
    % Set the groups from which you want to extract module specific summaries
        % only 4mo 
            cmask_4mo = [1,1; 0,0; 0,0];
    % set output file root for any outputs that are written out
        outroot11 = '11_CTRL_4MO_FEMALE_ref_ctrl_4mo_M';
    % set the reference partition
        ci_ctrl_4mo_M = ctrl_tier1.mrccPartition{2,2};
    % run the script
        fcn_summary_modules(ci_ctrl_4mo_M,ctrl_tier1,cmask_4mo,outroot11)        

% (12) CTRL 6mo FEMALE (REF = CTRL 6mo MALE)
    % Set the groups from which you want to extract module specific summaries
        % only 6mo 
            cmask_6mo = [0,0; 1,1; 0,0];
    % set output file root for any outputs that are written out
        outroot12 = '12_CTRL_6MO_FEMALE_ref_ctrl_6mo_M';
    % set the reference partition
        ci_ctrl_6mo_M = ctrl_tier1.mrccPartition{3,2};
    % run the script
        fcn_summary_modules(ci_ctrl_6mo_M,ctrl_tier1,cmask_6mo,outroot12) 
        
% (13) CTRL 12mo FEMALE (REF = CTRL 12mo MALE)
    % Set the groups from which you want to extract module specific summaries
        % only 12mo 
            cmask_12mo = [0,0; 0,0; 1,1];
    % set output file root for any outputs that are written out
        outroot13 = '13_CTRL_12MO_FEMALE_ref_ctrl_12mo_M';
    % set the reference partition
        ci_ctrl_12mo_M = ctrl_tier1.mrccPartition{4,2};
    % run the script
        fcn_summary_modules(ci_ctrl_12mo_M,ctrl_tier1,cmask_12mo,outroot13) 
        
% XFAD
% (14) XFAD 4mo FEMALE (REF = XFAD 4mo MALE)
    % Set the groups from which you want to extract module specific summaries
        % only 4mo 
            cmask_4mo = [1,1; 0,0; 0,0];
    % set output file root for any outputs that are written out
        outroot14 = '14_XFAD_4MO_FEMALE_ref_xfad_4mo_M';
    % set the reference partition
        ci_xfad_4mo_M = xfad_tier1.mrccPartition{2,2};
    % run the script
        fcn_summary_modules(ci_xfad_4mo_M,xfad_tier1,cmask_4mo,outroot14)        

% (15) XFAD 6mo FEMALE (REF = XFAD 6mo MALE)
    % Set the groups from which you want to extract module specific summaries
        % only 6mo 
            cmask_6mo = [0,0; 1,1; 0,0];
    % set output file root for any outputs that are written out
        outroot15 = '15_XFAD_6MO_FEMALE_ref_xfad_6mo_M';
    % set the reference partition
        ci_xfad_6mo_M = xfad_tier1.mrccPartition{3,2};
    % run the script
        fcn_summary_modules(ci_xfad_6mo_M,xfad_tier1,cmask_6mo,outroot15) 
        
% (16) XFAD 12mo FEMALE (REF = XFAD 12mo MALE)
    % Set the groups from which you want to extract module specific summaries
        % only 12mo 
            cmask_12mo = [0,0; 0,0; 1,1];
    % set output file root for any outputs that are written out
        outroot16 = '16_XFAD_12MO_FEMALE_ref_xfad_12mo_M';
    % set the reference partition
        ci_xfad_12mo_M = xfad_tier1.mrccPartition{4,2};
    % run the script
        fcn_summary_modules(ci_xfad_12mo_M,xfad_tier1,cmask_12mo,outroot16)
  
 
%% **************************** RS- HAVEN'T RUN PAST THIS YET... *************************************


%% tier 2 compare network metrics

     % Compare Males at each age (Set by the cmask input) between control
     % and xfad at p-value thresholds set in tier1.
     
     cmask_all = [1,1;1,1;1,1];

     % Note: Network metric comparisons are only done on thresholded networks. 
covariance_analysis_tier2a(ctrl_tier1,xfad_tier1,cmask_all,0.05)

 


% Compare within module network metrics (given a reference partition)
    
    
    %%         ci_xfad_12mo_M = xfad_tier1.mrccPartition{4,2};

     ci_ctrl_4mo_M = ctrl_tier1.mrccPartition{2,2};
    
covariance_analysis_tier2b(ci_ctrl_4mo_M,'CTRL_4mM_ref',ctrl_tier1,xfad_tier1,cmask_all,bluered_cmap,pval)












