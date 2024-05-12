%% Fig1J RdLight1 Saline vs Reserpine

% Choose and load the appropriate files
clear all;
[FPfiles,FPpath] = uigetfile('*.mat','Select FP Files to Analyze','MultiSelect','On');
if isempty(FPfiles)
    fprintf('No Photometry File Selected!\n');
else
    if (~iscell(FPfiles)); FPfiles = {FPfiles}; end
    load(fullfile(FPpath,FPfiles{1})); fprintf(['Processing File: ',FPfiles{1},'....\n\n']);  % Load each file
end

nTraces = 20;          % number of randomly selected traces to plot
% Plotting 4ms, so just include the following trials that correspond to 4ms conditions across all mice
idx = [31:45,31+60:45+60,31+120:45+120,31+180:45+180]';
if length(idx) > nTraces
    idx = datasample(idx,nTraces,'Replace',false); % choose N random trials
end

data.final(1).allFP_baseline_concat = [];
for x = 1:16
    for y = 1:15
        baseline_FP = data.final(x).allFP_dFF(:,y);       % Extracting each trial
        baseline_FP(20001:20111) = NaN;                   % stimulus mask
        baseline_FP = baseline_FP - nanmean(baseline_FP(1:20000)); % baselining
        baseline_FP = baseline_FP*100; % Convert decimal points into percentage
        data.final(1).allFP_baseline_concat = [data.final(1).allFP_baseline_concat,baseline_FP];
    end
end
% Smoothing/downsampling factor: 100
data.final(1).allFP_baseline_concat_downsample = data.final(1).allFP_baseline_concat(1:100:end,:);
time_downsample = 1:100:length(data.acq(1).time);

clr = {'m','r'}; % Used for bilateral RdLight (magenta- NAc, red- DLS)
nan_index = 201:202;
fig = figure; fig.Position(3) = 1000; clearvars sp
y = 1;
hold on;
plot(time_downsample/data.Fs, data.final(1).allFP_baseline_concat_downsample(:,idx), 'Color', [0.3 0.3 0.3],'LineWidth',0.25);       % For bilateral RdLight experiments

shadederrbar([time_downsample(1:200)/data.Fs]', nanmean(data.final(1).allFP_baseline_concat_downsample(1:200,idx),2), SEM(data.final(1).allFP_baseline_concat_downsample(1:200,idx),2), clr{y});    % For bilateral RdLight experiments
shadederrbar([time_downsample(203:end)/data.Fs]', nanmean(data.final(1).allFP_baseline_concat_downsample(203:end,idx),2), SEM(data.final(1).allFP_baseline_concat_downsample(203:end,idx),2), clr{y});    % For bilateral RdLight experiments

xlabel('time to stimulus (s)'); xlim([1 4]);
ylabel('fluorescence (%dF/F)');
title(sprintf('example (%d traces)',nTraces));

movegui(gcf,'center')