function [data] = JT_processFP_clipthenprocess(data,params)
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

% 230816 UPDATE: My updated script allows additional parameters for masking out
% the stimulation duration, such that it gets rid of the ripple-effect
% when the raw data is filtered/interpolated.

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
    % Make sure to average all of the corresponding hemisphere data before
    % any processing/filtering is applied.

    for x = 1:nFP
        rawFP = data.acq(n).FP{x};

        % 230815 UPDATE: Masks out stimulation duration:
        %rawFP(params.mask_win(1):params.mask_win(2),:)=NaN;
        %rawFP(params.mask_win(1):params.mask_win(2),:)=[];

    end
    
        nbFP = filterFP(rawFP,rawFs,lpCut,filtOrder,'lowpass');
        nbFP = downsampleTLab(nbFP,dsRate,dsType);
        [FP,baseline] = baselineFP(nbFP,interpType,fitType,basePrc,winSize,winOv,Fs);
        data.final(n).FP{x} = FP;
        data.final(n).nbFP{x} = nbFP;
        data.final(n).FPbaseline{x} = baseline;
        masking_vec=NaN(length(params.mask_win(1)/params.dsRate:params.mask_win(2)/params.dsRate),1);
        %data_test=[data.final(n).FP{x,1}(1:params.mask_win(1)-1); masking_vec; data.final(n).FP{x,1}(params.mask_win(2)+1:end)];
        data.final(n).masked_FP{x}=[data.final(n).FP{x,1}(1:round(params.mask_win(1)/params.dsRate)); masking_vec; data.final(n).FP{x,1}(round(params.mask_win(2)/params.dsRate):end)];
        data.final(n).masked_nbFP{x}=[data.final(n).nbFP{x,1}(1:round(params.mask_win(1)/params.dsRate)); masking_vec; data.final(n).nbFP{x,1}(round(params.mask_win(2)/params.dsRate):end)];
        data.final(n).masked_FPbaseline{x}=[data.final(n).FPbaseline{x,1}(1:round(params.mask_win(1)/params.dsRate)); masking_vec; data.final(n).FPbaseline{x,1}(round(params.mask_win(2)/params.dsRate):end)];
    end
    timeVec = [1:L]/Fs;
    data.final(n).time = timeVec';
end