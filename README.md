# CovNet
A workflow for analysis of covariance networks from brain imaging data.

Written in Matlab 2021a.
A description of outputs can be found in the <Covariance Analysis Workflow.docx> document, along with relevant refences.
 
Start by creating a project directory and downloading the CovNet repository into it. 
 Copy the covariance_analysis_workflow.m file out of CovNet and into the project directrory.
  Follow the comments in the workflow script to format your data for input into the script, then run the workflow script with the desired options.
 
 The code writes out figures in .pdf format, containing visualization of networks, summary measures, and comparisons. 

# Dependencies
 - Brain Connectivity Toolbox https://sites.google.com/site/bctnet/
 - GenLouvain Modularity Package https://github.com/GenLouvain/GenLouvain
 - Hierarchical Consensus Package https://github.com/LJeub/HierarchicalConsensus
 
These are currently included in this repository for ease of continued development, but will likely become external dependencies in later complete/stable versions.
