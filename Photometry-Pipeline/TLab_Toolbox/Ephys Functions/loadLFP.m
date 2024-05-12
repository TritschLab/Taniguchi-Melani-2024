function lfp = loadLFP(nChan,chan2pull,acqFs,dsFs)
%loadLFP - Load LFPs from In Vivo Ephys Recordings
%
%   Usage:
%       fp = loadLFP(nChan,chan2pull,acqFs,dsFs)
%
%   Description: This function pulls a single channel recording from an In
%   Vivo Ephys recording and decimates the trace to allow for LFP analysis
%   
%   Input:
%       nChan - Number of channels in the dat file
%       chan2pull - The single channel to pull from a recording
%       acqFs - The raw acquisition sampling frequency
%       dsFs - The new downsampled frequency%
%
%   Output:
%       lfp - The decimated channel data
%
%
%   Author: Pratik Mistry, 2020
%

[datFile,datPath] = uigetfile('*.dat','Select Dat File to Extract Single Channel LFP');
datFF = fullfile(datPath,datFile);
datMap = getDatMemMap(datFF,nChan);
lfp_full = double(datMap.Data.dat(chan2pull,:));
lfp = decimateTLab(lfp_full,acqFs,acqFs/dsFs);
end