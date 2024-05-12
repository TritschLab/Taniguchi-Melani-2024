%% RdLight Paper Figure 1E

%% Light Artifact Across Time: NAc
load('N1609A_fig1d_L.mat')
concatenate = [data.final(3).allFP_dFF(:,2); data.final(3).allFP_dFF(:,3); data.final(3).allFP_dFF(:,4); data.final(3).allFP_dFF(:,5); data.final(3).allFP_dFF(:,6)];
% Need to multiply the baselined values by 100 to get percentage (instead of fractional)
concatenate = concatenate*100;
% omit the light stimuli itself
plot_concatenate = concatenate;
plot_concatenate([20001:20040,90001:90040,160001:160040,230001:230040,300001:300040]) = NaN;

% Downsample by a factor of x100 down to 100Hz
smoothing_window = 100; smoothing_overlay = 0;
initial_index = []; final_index = [];
% AVERAGING/SMOOTHING
smoothed_interval = [];
% For the last stim, use the method of manually determining the interval of which to look at. (Default 4000ms interval)
for zz = 1:smoothing_window-smoothing_overlay:length(plot_concatenate)-smoothing_window
    smoothed_interval(zz,1) = mean(plot_concatenate(zz+round(smoothing_window/2)-round(smoothing_window/2):zz+round(smoothing_window/2)+round(smoothing_window/2)));
end
smoothed_interval = [smoothed_interval(smoothed_interval~=0)];

figure(1);
plot((1:length(smoothed_interval))/smoothing_window,smoothed_interval)
ylim([-4 8])
xlim([0 35])
xlabel('seconds')
ylabel('%dF/F')
title('NAc 5 Pulses @ 9mW, 4ms duration (N1609A-fig1d-NAcOnly)')
box off


%% Light Artifact Across Time: DLS
load('N1609A_fig1d_R.mat')
concatenate = [data.final(3).allFP_dFF(:,2); data.final(3).allFP_dFF(:,3); data.final(3).allFP_dFF(:,4); data.final(3).allFP_dFF(:,5); data.final(3).allFP_dFF(:,6)];
% plot(concatenate)
% Need to multily the baselined values by 100 to get percentage (instead of fractional)
concatenate = concatenate*100;
% omit the light stimuli itself
plot_concatenate = concatenate;
plot_concatenate([20001:20040,90001:90040,160001:160040,230001:230040,300001:300040]) = NaN;

% Downsample by a factor of x100 down to 100Hz
smoothing_window = 100; smoothing_overlay = 0;
initial_index = []; final_index = [];
% AVERAGING/SMOOTHING
smoothed_interval = [];
% For the last stim, use the method of manually determining the interval of which to look at. (Default 4000ms interval)
for zz = 1:smoothing_window-smoothing_overlay:length(plot_concatenate)-smoothing_window
    smoothed_interval(zz,1) = mean(plot_concatenate(zz+round(smoothing_window/2)-round(smoothing_window/2):zz+round(smoothing_window/2)+round(smoothing_window/2)));
end
smoothed_interval = [smoothed_interval(smoothed_interval~=0)];

figure(2);
plot((1:length(smoothed_interval))/smoothing_window,smoothed_interval)
ylim([-4 8])
xlim([0 35])
xlabel('seconds')
ylabel('%dF/F')
title('DLS 5 Pulses @ 9mW, 4ms duration (N1609A-fig1d-NAcOnly)')
box off
