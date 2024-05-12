%% Fig1G_Fig2B_Fig2C

% Choose and load the appropriate files for analysis
clear all;
[FPfiles,FPpath] = uigetfile('*.mat','Select FP Files to Analyze','MultiSelect','On');

%% Fig1d-rDA3m mean response across time

% NOTE: Using JT065, Right-Hemi NAc as the example mouse (I only did n = 1 for this experiment)

% Number of conditions
process_idx = 4;

% Initialize Plot Colors: Using colors that Berke used in his plot
clr={'#351111','#620c0d','#b31f24','#ed2124'};
% 1ms: 351111/[53 17 17]
% 2ms: 620c0d/[98 12 13]
% 4ms: b31f24/[179 31 36]
% 6ms: ed2124/[237 33 
% 36]
for i=1:length(clr)
    str=clr{i};
    color(i,:) = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
end

% Clip out the stimuli itself x3.0 except for 1ms (x5)
pulse_width_vec = [0.001 0.002 0.004 0.006];
for y = 1:process_idx
    if (y == 1)
        data.final(y).baselineFP_dFF(data.onset_idx(1,y):data.offset_idx(1,y)+pulse_width_vec(y)*data.Fs*4.0) = NaN;
        data.final(y).SEM_FP_dFF(data.onset_idx(1,y):data.offset_idx(1,y)+pulse_width_vec(y)*data.Fs*4.0) = NaN;
        initial_index(y) = data.onset_idx(1,y)-1;
        final_index(y) = data.offset_idx(1,y)+pulse_width_vec(y)*data.Fs*4.0+1;
    else
        data.final(y).baselineFP_dFF(data.onset_idx(1,y):data.offset_idx(1,y)+pulse_width_vec(y)*data.Fs*2.0) = NaN;
        data.final(y).SEM_FP_dFF(data.onset_idx(1,y):data.offset_idx(1,y)+pulse_width_vec(y)*data.Fs*2.0) = NaN;
        initial_index(y) = data.onset_idx(1,y)-1;
        final_index(y) = data.offset_idx(1,y)+pulse_width_vec(y)*data.Fs*2.0+1;
    end

end

f = figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
for y = 1:process_idx
    % Plotting the mean dFF
    % Superimposing
    subplot(process_idx,1,y)
    shadederrbar(data.acq(1).time(1:initial_index(y)), data.final(y).baselineFP_dFF(1:initial_index(y)), data.final(y).SEM_FP_dFF(1:initial_index(y)), color(y,:)); % plot average across trials
    shadederrbar(data.acq(1).time(final_index(y):end), data.final(y).baselineFP_dFF(final_index(y):end), data.final(y).SEM_FP_dFF(final_index(y):end), color(y,:)); % plot average across trials
    %shadederrbar(data.acq(1).time(data.final(y).baselineFP_dFF~=0), data.final(y).baselineFP_dFF(data.final(y).baselineFP_dFF~=0), data.final(y).SEM_FP_dFF(data.final(y).SEM_FP_dFF~=0), color(y,:)); % plot average across trials

    if (y == 1)
        title('NAc rDA3m Fluorescence','FontSize',40);
        %title('DLS rDA3m Fluorescence','FontSize',40);
    end
    xlabel('Seconds','FontSize',30);
    ylabel('%dF/F','FontSize',30);
    xlim([1.75 3]);
    ylim([-4 3]);
    H=gca;
    H.LineWidth=5; %change to the desired value 
end
hold off;

% Save the figure
saveas(f,strcat('JT050 rDA3m Fig1d leftNAc 9mW Raw Mean Traces from n = 1','.fig')); % .fig
saveas(f,strcat('JT050 rDA3m Fig1d leftNAc 9mW Raw Mean Traces from n = 1','.png')); % .png
%saveas(f,strcat('JT065 rDA3m Fig1d DLS 9mW Raw Mean Traces from n = 1','.fig')); % .fig
%saveas(f,strcat('JT065 rDA3m Fig1d DLS 9mW Raw Mean Traces from n = 1','.png')); % .png

f = figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
for y = 1:process_idx
    % Plotting the mean dFF
    % Superimposing
    subplot(process_idx,1,y)
    shadederrbar(data.acq(1).time(1:initial_index(y)), data.final(y).baselineFP_dFF(1:initial_index(y)), data.final(y).SEM_FP_dFF(1:initial_index(y)), color(y,:)); % plot average across trials
    shadederrbar(data.acq(1).time(final_index(y):end), data.final(y).baselineFP_dFF(final_index(y):end), data.final(y).SEM_FP_dFF(final_index(y):end), color(y,:)); % plot average across trials
    %shadederrbar(data.acq(1).time(data.final(y).baselineFP_dFF~=0), data.final(y).baselineFP_dFF(data.final(y).baselineFP_dFF~=0), data.final(y).SEM_FP_dFF(data.final(y).SEM_FP_dFF~=0), color(y,:)); % plot average across trials


    if (y == 1)
        title('NAc rDA3m Fluorescence','FontSize',40);
        %title('DLS rDA3m Fluorescence','FontSize',40);
    end
    xlabel('Seconds','FontSize',30);
    ylabel('%dF/F','FontSize',30);
    xlim([1.75 3]);
    ylim([-4 3]);
    H=gca;
    H.LineWidth=5; %change to the desired value 
    set(gca,'xticklabel',[]);
    set(gca,'yticklabel',[]);
end
hold off;

% Save the figure
saveas(f,strcat('JT050 rDA3m Fig1d leftNAc 9mW Raw Mean Traces from n = 1-no-ticks','.fig')); % .fig
saveas(f,strcat('JT050 rDA3m Fig1d leftNAc 9mW Raw Mean Traces from n = 1-no-ticks','.png')); % .png
%saveas(f,strcat('JT065 rDA3m Fig1d DLS 9mW Raw Mean Traces from n = 1-no-ticks','.fig')); % .fig
%saveas(f,strcat('JT065 rDA3m Fig1d DLS 9mW Raw Mean Traces from n = 1-no-ticks','.png')); % .png

%% Fig1G

% NOTE: Using N1609B-9mW NAc/DLS as example

% Number of conditions
process_idx = 4;

% Initialize Plot Colors: Using colors that Berke used in his plot
clr={'#351111','#620c0d','#b31f24','#ed2124'};
% 1ms: 351111/[53 17 17]
% 2ms: 620c0d/[98 12 13]
% 4ms: b31f24/[179 31 36]
% 6ms: ed2124/[237 33 36]
for i=1:length(clr)
    str=clr{i};
    color(i,:) = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
end

% Clip out the stimuli itself x2.0 except for 1ms (x3)
pulse_width_vec = [0.001 0.002 0.004 0.006];
for y = 1:process_idx
    if (y == 1)
        data.final(y).baselineFP_dFF(data.onset_idx(1,y):data.offset_idx(1,y)+pulse_width_vec(y)*data.Fs*2.0) = NaN;
        data.final(y).SEM_FP_dFF(data.onset_idx(1,y):data.offset_idx(1,y)+pulse_width_vec(y)*data.Fs*2.0) = NaN;
    else
        data.final(y).baselineFP_dFF(data.onset_idx(1,y):data.offset_idx(1,y)+pulse_width_vec(y)*data.Fs*1.0) = NaN;
        data.final(y).SEM_FP_dFF(data.onset_idx(1,y):data.offset_idx(1,y)+pulse_width_vec(y)*data.Fs*1.0) = NaN;
    end

end

data.smoothing_window = 100;
data.smoothing_overlay = 0;

for y = 1:process_idx

        % AVERAGING/SMOOTHING
        % Default: 1ms window size. 0.5ms overlay.
        smoothed_interval = [];
        smoothed_interval_sem = [];
        smoothed_interval_time = [];
        % For the last stim, use the method of manually determining the interval of which to look at. (Default 4000ms interval)
        extracted_interval = data.final(y).baselineFP_dFF;
        extracted_interval_sem = data.final(y).SEM_FP_dFF;
        extracted_interval_time = data.acq(1).time;

        for zz = 1:data.smoothing_window-data.smoothing_overlay:length(extracted_interval)-data.smoothing_window
            smoothed_interval(zz,1) = mean(extracted_interval(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
            smoothed_interval_sem(zz,1) = mean(extracted_interval_sem(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
            smoothed_interval_time(zz,1) = mean(extracted_interval_time(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
        end

        % Align the Smoothed Interval to Real Intervals (i.e. the time from start of interval to time of 1st sampling) for precise indexing
        data.final(y).all_smoothed_interval_raw(:,1) = smoothed_interval;
        data.final(y).all_smoothed_interval_sem_raw(:,1) = smoothed_interval_sem;
        data.final(y).all_smoothed_interval_time_raw(:,1) = smoothed_interval_time;
        smoothed_interval_work = [smoothed_interval(smoothed_interval~=0)];
        smoothed_interval_work_sem = [smoothed_interval_sem(smoothed_interval_sem~=0)];
        smoothed_interval_work_time = [smoothed_interval_time(smoothed_interval_time~=0)];
        data.final(y).all_smoothed_interval_idx_raw(:,1) = [find(smoothed_interval~=0)+data.smoothing_window/2];
        data.final(y).all_smoothed_interval_idx_sem_raw(:,1) = [find(smoothed_interval_sem~=0)+data.smoothing_window/2];
        data.final(y).all_smoothed_interval_idx_time_raw(:,1) = [find(smoothed_interval_time~=0)+data.smoothing_window/2];
        % Remove all values that are the stimulus itself (Approx Stimuli Interval * 3)
        % Remove all values related with the stimulus fluctuations (i.e. very high numbers and low numbers)
        smoothed_interval_work(smoothed_interval_work>20) = 0;
        smoothed_interval_work_sem(smoothed_interval_work_sem>20) = 0;
        smoothed_interval_work(smoothed_interval_work<-5) = 0;
        smoothed_interval_work_sem(smoothed_interval_work_sem<-5) = 0;
        % Remove all zeros (i.e. non-sampled indicies from the vector)
        smoothed_interval_work(smoothed_interval_work == 0) = NaN;
        smoothed_interval_work_sem(smoothed_interval_work_sem == 0) =NaN;
        data.final(y).all_smoothed_interval_processed(:,1) = smoothed_interval_work;
        data.final(y).all_smoothed_interval_sem_processed(:,1) = smoothed_interval_work_sem;
        data.final(y).all_smoothed_interval_time_processed(:,1) = smoothed_interval_work_time;

        initial_index(y) = find(isnan(smoothed_interval_work),1,'first')-1;
        final_index(y) = find(isnan(smoothed_interval_work),1,'last')+1;
end

f = figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
for y = 1:process_idx
    % Plotting the mean dFF
    % Superimposing
    shadederrbar(data.final(y).all_smoothed_interval_time_processed(1:initial_index(y)), data.final(y).all_smoothed_interval_processed(1:initial_index(y)), data.final(y).all_smoothed_interval_sem_processed(1:initial_index(y)), color(y,:)); % plot average across trials
    shadederrbar(data.final(y).all_smoothed_interval_time_processed(final_index(y):end), data.final(y).all_smoothed_interval_processed(final_index(y):end), data.final(y).all_smoothed_interval_sem_processed(final_index(y):end), color(y,:)); % plot average across trials
end
hold off;
xlim([1 6]);
ylim([-5 10]);
H=gca;
H.LineWidth=5; %change to the desired value   

% Save the figure with tick labels
saveas(f,strcat('Fig1d-RdLight-N1609B-2mW-DLS-Mean-Traces-from-N=1','.fig')); % .fig
saveas(f,strcat('Fig1d-RdLight-N1609B-2mW-DLS-Mean-Traces-from-N=1','.png')); % .png
%saveas(f,strcat('Fig1d-RdLight-N1609B-DLS-Mean-Traces-from-N=1','.fig')); % .fig
%saveas(f,strcat('Fig1d-RdLight-N1609B-DLS-Mean-Traces-from-N=1','.png')); % .png

set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);

% Save the figure without tick labels
saveas(f,strcat('Fig1d-RdLight-N1609B-2mW-DLS-Mean-Traces-from-N=1-without-ticks','.fig')); % .fig
saveas(f,strcat('Fig1d-RdLight-N1609B-2mW-DLS-Mean-Traces-from-N=1-without-ticks','.png')); % .png
%saveas(f,strcat('Fig1d-RdLight-N1609B-DLS-Mean-Traces-from-N=1-without-ticks','.fig')); % .fig
%saveas(f,strcat('Fig1d-RdLight-N1609B-DLS-Mean-Traces-from-N=1-without-ticks','.png')); % .png

%% Fig2B

% NOTE: Using N1609A-9mW NAc as example

% Number of conditions
process_idx = 5;

% Initialize Plot Colors: Using colors that Berke used in his plot
clr={'#030404','#351111','#620c0d','#b31f24','#ed2124'};
% 1p: 030404/[3 4 4]
% 4p: 351111/[53 17 17]
% 16p: 620c0d/[98 12 13]
% 32p: b31f24/[179 31 36]
% 64p: ed2124/[237 33 36]
for i=1:length(clr)
    str=clr{i};
    color(i,:) = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
end

% Clip out the stimuli itself x2.0
for y = 1:process_idx
    if (y == 1)
        data.final(y).baselineFP_dFF(data.onset_idx(1,y):data.offset_idx(1,y)+data.pulse_width(1)*data.Fs*1.0) = NaN;
        data.final(y).SEM_FP_dFF(data.onset_idx(1,y):data.offset_idx(1,y)+data.pulse_width(1)*data.Fs*1.0) = NaN;
    else
        for z = 1 : nnz(~isnan(data.onset_idx(:,y)))
            data.final(y).baselineFP_dFF(data.onset_idx(z,y):data.offset_idx(z,y)+data.pulse_width(1)*data.Fs*1.0) = NaN;
            data.final(y).SEM_FP_dFF(data.onset_idx(z,y):data.offset_idx(z,y)+data.pulse_width(1)*data.Fs*1.0) = NaN;
        end
    end
end

data.smoothing_window = 50;
data.smoothing_overlay = 0;

for y = 1:process_idx

        % AVERAGING/SMOOTHING
        % Default: 1ms window size. 0.5ms overlay.
        smoothed_interval = [];
        smoothed_interval_sem = [];
        smoothed_interval_time = [];
        % For the last stim, use the method of manually determining the interval of which to look at. (Default 4000ms interval)
        extracted_interval = data.final(y).baselineFP_dFF;
        extracted_interval_sem = data.final(y).SEM_FP_dFF;
        extracted_interval_time = data.acq(1).time;

        for zz = 1:data.smoothing_window-data.smoothing_overlay:length(extracted_interval)-data.smoothing_window
            smoothed_interval(zz,1) = mean(extracted_interval(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
            smoothed_interval_sem(zz,1) = mean(extracted_interval_sem(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
            smoothed_interval_time(zz,1) = mean(extracted_interval_time(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
        end

        % Align the Smoothed Interval to Real Intervals (i.e. the time from start of interval to time of 1st sampling) for precise indexing
        data.final(y).all_smoothed_interval_raw(:,1) = smoothed_interval;
        data.final(y).all_smoothed_interval_sem_raw(:,1) = smoothed_interval_sem;
        data.final(y).all_smoothed_interval_time_raw(:,1) = smoothed_interval_time;
        smoothed_interval_work = [smoothed_interval(smoothed_interval~=0)];
        smoothed_interval_work_sem = [smoothed_interval_sem(smoothed_interval_sem~=0)];
        smoothed_interval_work_time = [smoothed_interval_time(smoothed_interval_time~=0)];
        data.final(y).all_smoothed_interval_idx_raw(:,1) = [find(smoothed_interval~=0)+data.smoothing_window/2];
        data.final(y).all_smoothed_interval_idx_sem_raw(:,1) = [find(smoothed_interval_sem~=0)+data.smoothing_window/2];
        data.final(y).all_smoothed_interval_idx_time_raw(:,1) = [find(smoothed_interval_time~=0)+data.smoothing_window/2];
        % Remove all values that are the stimulus itself (Approx Stimuli Interval * 3)
        % Remove all values related with the stimulus fluctuations (i.e. very high numbers and low numbers)
        smoothed_interval_work(smoothed_interval_work>20) = 0;
        smoothed_interval_work_sem(smoothed_interval_work_sem>20) = 0;
        smoothed_interval_work(smoothed_interval_work<-5) = 0;
        smoothed_interval_work_sem(smoothed_interval_work_sem<-5) = 0;
        % Remove all zeros (i.e. non-sampled indicies from the vector)
        smoothed_interval_work(smoothed_interval_work == 0) = [];
        smoothed_interval_work_sem(smoothed_interval_work_sem == 0) = [];
        data.final(y).all_smoothed_interval_processed(:,1) = smoothed_interval_work;
        data.final(y).all_smoothed_interval_sem_processed(:,1) = smoothed_interval_work_sem;
        data.final(y).all_smoothed_interval_time_processed(:,1) = smoothed_interval_work_time;

        total_index = find(isnan(smoothed_interval_work));
        mask_counter = 0;
        mask_loop = 0;
        for zz = 1:length(total_index)
            if (mask_loop == 0)
                if (mask_counter == 0)
                     initial_index(y,1) = total_index(zz);
                else
                     initial_index(y,nnz(initial_index(y,:))+1) = total_index(zz);
                end
                mask_loop = 1;
                mask_counter = mask_counter + 1;
            else
                if (zz ~= length(total_index))
                  if (total_index(zz+1)-total_index(zz)>1)
                      if (nnz(initial_index(y,:)) == 1)
                          final_index(y,1) = total_index(zz);
                      else
                          final_index(y,nnz(initial_index(y,:))) = total_index(zz);
                      end
                      mask_loop = 0;
                  end
                else
                    final_index(y,nnz(initial_index(y,:))) = total_index(zz);
                end
            end   
        end

end

f = figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
for y = 1:process_idx
    % Plotting the mean dFF
    % Superimposing
    for z = 1:nnz(initial_index(y,:))
      shadederrbar(data.final(y).all_smoothed_interval_time_processed(1:initial_index(y,z)-1), data.final(y).all_smoothed_interval_processed(1:initial_index(y,z)-1), data.final(y).all_smoothed_interval_sem_processed(1:initial_index(y,z)-1), color(y,:)); % plot average across trials
      shadederrbar(data.final(y).all_smoothed_interval_time_processed(final_index(y,z)+1:end), data.final(y).all_smoothed_interval_processed(final_index(y,z)+1:end), data.final(y).all_smoothed_interval_sem_processed(final_index(y,z)+1:end), color(y,:)); % plot average across trials
    end
end
hold off;
xlim([0 10]);
ylim([-5 25]);
H=gca;
H.LineWidth=5; %change to the desired value   

% Save the figure with tick labels
%saveas(f,strcat('Fig1c-RdLight-N1609A-NAc-Mean-Traces-from-N=1','.fig')); % .fig
%saveas(f,strcat('Fig1c-RdLight-N1609A-NAc-Mean-Traces-from-N=1','.png')); % .png
saveas(f,strcat('Fig1c-RdLight-N1609A-2mW-DLS-Mean-Traces-from-N=1','.fig')); % .fig
saveas(f,strcat('Fig1c-RdLight-N1609A-2mW-DLS-Mean-Traces-from-N=1','.png')); % .png

set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);

% Save the figure without tick labels
%saveas(f,strcat('Fig1c-RdLight-N1609A-NAc-Mean-Traces-from-N=1-without-ticks','.fig')); % .fig
%saveas(f,strcat('Fig1c-RdLight-N1609A-NAc-Mean-Traces-from-N=1-without-ticks','.png')); % .png
saveas(f,strcat('Fig1c-RdLight-N1609A-2mW-DLS-Mean-Traces-from-N=1-without-ticks','.fig')); % .fig
saveas(f,strcat('Fig1c-RdLight-N1609A-2mW-DLS-Mean-Traces-from-N=1-without-ticks','.png')); % .png

%% Fig2C

% NOTE: Using JT070-9mW NAc as example

% Number of conditions
process_idx = 3;

% Initialize Plot Colors: Using colors that Berke used in his plot
clr={'#ed2124','#b31f24','#620c0d'};
% 2.5s: 620c0d/[3 4 4]
% 5s: b31f24/[53 17 17]
% 10s: ed2124/[98 12 13]
for i=1:length(clr)
    str=clr{i};
    color(i,:) = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
end

% Clip out the stimuli itself x2.5
for y = 1:process_idx
        for z = 1 : nnz(~isnan(data.onset_idx(:,y)))
            data.final(y).baselineFP_dFF(data.onset_idx(z,y):data.offset_idx(z,y)+0.004*data.Fs*1.5) = NaN;
            data.final(y).SEM_FP_dFF(data.onset_idx(z,y):data.offset_idx(z,y)+0.004*data.Fs*1.5) = NaN;
        end
end

data.smoothing_window = 50;
data.smoothing_overlay = 0;

for y = 1:process_idx

        % AVERAGING/SMOOTHING
        % Default: 1ms window size. 0.5ms overlay.
        smoothed_interval = [];
        smoothed_interval_sem = [];
        smoothed_interval_time = [];
        % For the last stim, use the method of manually determining the interval of which to look at. (Default 4000ms interval)
        extracted_interval = data.final(y).baselineFP_dFF;
        extracted_interval_sem = data.final(y).SEM_FP_dFF;
        extracted_interval_time = data.acq(1).time;

        for zz = 1:data.smoothing_window-data.smoothing_overlay:length(extracted_interval)-data.smoothing_window
            smoothed_interval(zz,1) = mean(extracted_interval(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
            smoothed_interval_sem(zz,1) = mean(extracted_interval_sem(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
            smoothed_interval_time(zz,1) = mean(extracted_interval_time(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
        end

        % Align the Smoothed Interval to Real Intervals (i.e. the time from start of interval to time of 1st sampling) for precise indexing
        data.final(y).all_smoothed_interval_raw(:,1) = smoothed_interval;
        data.final(y).all_smoothed_interval_sem_raw(:,1) = smoothed_interval_sem;
        data.final(y).all_smoothed_interval_time_raw(:,1) = smoothed_interval_time;
        smoothed_interval_work = [smoothed_interval(smoothed_interval~=0)];
        smoothed_interval_work_sem = [smoothed_interval_sem(smoothed_interval_sem~=0)];
        smoothed_interval_work_time = [smoothed_interval_time(smoothed_interval_time~=0)];
        data.final(y).all_smoothed_interval_idx_raw(:,1) = [find(smoothed_interval~=0)+data.smoothing_window/2];
        data.final(y).all_smoothed_interval_idx_sem_raw(:,1) = [find(smoothed_interval_sem~=0)+data.smoothing_window/2];
        data.final(y).all_smoothed_interval_idx_time_raw(:,1) = [find(smoothed_interval_time~=0)+data.smoothing_window/2];
        % Remove all values that are the stimulus itself (Approx Stimuli Interval * 3)
        % Remove all values related with the stimulus fluctuations (i.e. very high numbers and low numbers)
        smoothed_interval_work(smoothed_interval_work>20) = 0;
        smoothed_interval_work_sem(smoothed_interval_work_sem>20) = 0;
        smoothed_interval_work(smoothed_interval_work<-5) = 0;
        smoothed_interval_work_sem(smoothed_interval_work_sem<-5) = 0;
        % Remove all zeros (i.e. non-sampled indicies from the vector)
        smoothed_interval_work(smoothed_interval_work == 0) = NaN;
        smoothed_interval_work_sem(smoothed_interval_work_sem == 0) = NaN;
        data.final(y).all_smoothed_interval_processed(:,1) = smoothed_interval_work;
        data.final(y).all_smoothed_interval_sem_processed(:,1) = smoothed_interval_work_sem;
        data.final(y).all_smoothed_interval_time_processed(:,1) = smoothed_interval_work_time;

        total_index = find(isnan(smoothed_interval_work));
        mask_counter = 0;
        mask_loop = 0;
        for zz = 1:length(total_index)
            if (mask_loop == 0)
                if (mask_counter == 0)
                     initial_index(y,1) = total_index(zz);
                else
                     initial_index(y,nnz(initial_index(y,:))+1) = total_index(zz);
                end
                mask_loop = 1;
                mask_counter = mask_counter + 1;
            else
                if (zz ~= length(total_index))
                  if (total_index(zz+1)-total_index(zz)>1)
                      if (nnz(initial_index(y,:)) == 1)
                          final_index(y,1) = total_index(zz);
                      else
                          final_index(y,nnz(initial_index(y,:))) = total_index(zz);
                      end
                      mask_loop = 0;
                  end
                else
                    final_index(y,nnz(initial_index(y,:))) = total_index(zz);
                end
            end   
        end
end

f = figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
for y = 1:process_idx
    % Plotting the mean dFF
    % Superimposing
    for z = 1:nnz(initial_index(y,:))
      shadederrbar(data.final(y).all_smoothed_interval_time_processed(1:initial_index(y,z)-1), data.final(y).all_smoothed_interval_processed(1:initial_index(y,z)-1), data.final(y).all_smoothed_interval_sem_processed(1:initial_index(y,z)-1), color(y,:)); % plot average across trials
          if (z ~= nnz(initial_index(y,:)))
                shadederrbar(data.final(y).all_smoothed_interval_time_processed(final_index(y,z)+1:initial_index(y,z+1)-1), data.final(y).all_smoothed_interval_processed(final_index(y,z)+1:initial_index(y,z+1)-1), data.final(y).all_smoothed_interval_sem_processed(final_index(y,z)+1:initial_index(y,z+1)-1), color(y,:)); % plot average across trials
          else
                shadederrbar(data.final(y).all_smoothed_interval_time_processed(final_index(y,z)+1:end), data.final(y).all_smoothed_interval_processed(final_index(y,z)+1:end), data.final(y).all_smoothed_interval_sem_processed(final_index(y,z)+1:end), color(y,:)); % plot average across trials
          end
    end
end
hold off;
xlim([0 14]);
ylim([-2 10]);
H=gca;
H.LineWidth=5; %change to the desired value   

% Save the figure with tick labels
saveas(f,strcat('Fig1f-RdLight-JT070-2mW-DLS-Mean-Traces-from-N=1','.fig')); % .fig
saveas(f,strcat('Fig1f-RdLight-JT070-2mW-DLS-Mean-Traces-from-N=1','.png')); % .png

set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);

% Save the figure without tick labels
saveas(f,strcat('Fig1f-RdLight-JT070-2mW-DLS-Mean-Traces-from-N=1-without-ticks','.fig')); % .fig
saveas(f,strcat('Fig1f-RdLight-JT070-2mW-DLS-Mean-Traces-from-N=1-without-ticks','.png')); % .png
