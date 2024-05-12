function [data] = JT_processFP(data,params)
%Process Fiber Photometry
%
%   [data] = processFP(data,params)
%
%   Description: This function is designed to process fiber photometry data
%   for the lab. The function performs demodulation (if selected),
%   filtering, baselining, and downsampling for all photometry traces in
%   the recording. The parameters for the analysis are found in the params
%   structure, which is created from a user-created scripted based on the
%   processParam.m file.
%
%   Input:
%   - data - A data structure specific to the Tritsch Lab. Created using
%   the convertH5_FP script
%   - params - A structure created from a variant of the processParams
%   script
%
%   Output:
%   - data - Updated data structure containing processed data
%
%   Author: Pratik Mistry 2019

nAcq = length(data.acq);
lpCut = params.FP.lpCut; filtOrder = params.FP.filtOrder;

dsRate = params.dsRate;
dsType = params.dsType;

interpType = params.FP.interpType;
fitType = params.FP.fitType; winSize = params.FP.winSize;
winOv = params.FP.winOv;
basePrc = params.FP.basePrc;

rawFs = data.gen.acqFs;
Fs = rawFs/dsRate;
data.gen.Fs = Fs;

for n = 1:nAcq
    L = size(data.acq(n).time,1);
    nFP = data.acq(n).nFPchan;
    FPnames = data.acq(n).FPnames;
    data.final(n).FPnames = FPnames;
    data.final(n).nFPchan = nFP;
    data.final(n).FP = cell(nFP,1);
    data.final(n).nbFP = cell(nFP,1);
    data.final(n).FPbaseline = cell(nFP,1);
    L = length(1:dsRate:L);
    for x = 1:nFP
        try
            rawFP = data.acq(n).FP{x};
        catch
            rawFP = data.acq(n).FP;
        end
        % 240220 NOTE: Haixin filters the data, then caluclates the baseline (10 prc of bottom), and then does the downsampling
        % Different from Anya/Pratik's protocol that does the downsampling, then the baselining (which will choose fewer data points for baselining itself)
        nbFP = filterFP(rawFP,rawFs,lpCut,filtOrder,'lowpass');
        nbFP = downsampleTLab(nbFP,dsRate,dsType);
        [FP,baseline] = baselineFP(nbFP,interpType,fitType,basePrc,winSize,winOv,Fs);
        data.final(n).FP{x} = FP;
        data.final(n).nbFP{x} = nbFP;
        data.final(n).FPbaseline{x} = baseline;
    end
    timeVec = [1:L]/Fs;
    data.final(n).time = timeVec';
end