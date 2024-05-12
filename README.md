The MATLAB functions and script in this repository are developed to process, compute, and visualize data acquired during fiber photometry experiments presented in Taniguchi et al. 2024.

## Dependencies:

Make sure you add to your path the following folders:

1. ~/Photometry-Pipeline

2. ~/Behavioral-Pipeline

Contains essential functions for processing and plotting data that is aligned to the reward onset.

For more information, please look at these following repositories

Tritsch Lab Toolbox: https://github.com/TritschLab/TLab_Toolbox

Tritsch Lab Photometry Pipeline: https://github.com/pratikmistry96/Photometry-Pipeline

## Steps:

For reward-aligned data:

1. Convert h5 files to .mat files
   
       >> convertH5
2. Edit Parameters
   
       >> edit processParams_mod

3. Process Data
   
       >> processData
4. Plot the figures
   
       >> Fig1D_alignment_to_rewards
	   
For stimulation-aligned data:

1. Process the data
   
       >> JT_RdLight_analysis_combined_peakdip
   
2. Run the corresponding plotting code located in the home directory


