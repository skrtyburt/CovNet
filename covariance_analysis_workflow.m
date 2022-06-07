
% Prepare data:
% example:
% control animals at 3 ages and two sexes.
% Currently organized as animal by region.

% MAKE SURE THAT DATA THAT GOES INTO THE CELL ARE 
% REGION (ROWS) BY ANIMAL (COLUMN).
%%
% Example of how cell structure was buld:
%{

% load data that was imported from excel 
load('CTRL.mat')

% Organize data into a cell and save it:
%initialize cell
CTRL = cell.empty;
% group title
CTRL{1,1}='Control';
% row labels (here by age)
CTRL{2,1}='4m';
CTRL{3,1}='6m';
CTRL{4,1}='12m';
% column labels (here by sex)
CTRL{1,2}='M';
CTRL{1,3}='F';
% copy data into the corresponding {row,column} indices
% HERE I AM TRANSPOSING THE DATA (') SO THAT ROWS ARE REGIONS
CTRL{2,2}=CTRL_4M';  CTRL{2,3}=CTRL_4F';
CTRL{3,2}=CTRL_6M';  CTRL{3,3}=CTRL_6F';
CTRL{4,2}=CTRL_12M'; CTRL{4,3}=CTRL_12F';
% saving to an exisisting file (keeping the raw data and cellData in the
% same file).
save('CTRL.mat','CTRL','-append')

%}
%%
% add code package to path
addpath(genpath([pwd '/functions']))
load colormaps.mat

% load ROI labels
load('roi_labels.mat')

% define a single or a range of p-values for thresholding correlation
% significance within each network
pval = [.005];

% load and run control data
load('CTRL.mat')
covariance_analysis_tier1(CTRL,roi_labels,bluered_cmap,pval)

% load and run mutant data
load('XFAD.mat')
covariance_analysis_tier1(XFAD,roi_labels,bluered_cmap,pval)

%%
% set your reference modularity partition
    % example 4m M Control
    ci = ctrl_tier1.mrccPartition{2,2};

% set the groups from which you want to extract module specific summaries
    % example only the males
    cmask = [1,0; 1,0; 1,0];
    % this variable will look like this:
    % 1 0
    % 1 0
    % 1 0
    % it should be the same size as number of data cells (in our example
    % its 3 ages on rows and sex as two columns).
%%
% load tier1 output into a structure
ctrl_tier1 = load('covariance_out_tier1_Control.mat');
% set output file root for any outputs that are written out;
outroot1 = 'ctrl_4M_M_ref_ctrl';
% run the script
% function for extracting mean SUVR for specified input modules
% means and stdev for each module are saved to a text file
fcn_summary_modules(ci,ctrl_tier1,cmask,outroot1)
%%

% repeat for xfad using the same reference partition and cmask
xfad_tier1 = load('covariance_out_tier1_xfad.mat');
% if you want to change reference to xfad 4m M
ci2 = xfad_tier1.mrccPartition{2,2};
% label for the text file
outroot2 = 'xfad_4M_M_ref_xfad';
% run the script
fcn_summary_modules(ci2,xfad_tier1,cmask,outroot2)

%% tier 2 compare network metrics


ci3 = ctrl_tier1.mrccPartition{2,2};
cmask2 = [1,0;0,0;0,0];
covariance_analysis_tier2(ci3,'ctrl_4m_M',ctrl_tier1,xfad_tier1,cmask2)












