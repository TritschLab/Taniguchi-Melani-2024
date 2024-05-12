%% RdLight_analysis_combined

% Choose and load the appropriate files for analysis
clear all;
[FPfiles,FPpath] = uigetfile('*.mat','Select FP Files to Analyze','MultiSelect','On');

if isempty(FPfiles)
fprintf('No Photometry File Selected!\n');
else
% If any photometry (.mat) files are selected, then proceed the analysis
% Variable 'FP files' will store all of the filenames needed
if (~iscell(FPfiles))
    FPfiles = {FPfiles};
end
nFiles = length(FPfiles);   % Number of files currently on hand

% For each .mat file, proceed to concatenate all sweeps into one matrix
for x = 1:nFiles

    % Load each file
    load(fullfile(FPpath,FPfiles{x}));
    fprintf(['Processing File: ',FPfiles{x},'....\n\n']);

    % Parameters
    if (strfind(FPfiles{1,1},'L'))
        data.location = "NAc-medial-core";
    elseif strfind(FPfiles{1,1},'R')
        data.location = "DLS";
    elseif strfind(FPfiles{1,1},'H')
        data.location = "HPC";
    end
    
    % Collect and group all data per hemisphere (i.e. region)
    % Which channel corresponds to the RdLight (i.e. target) and stimulus (i.e. blue light) channels of interest
    data.target_FPchannel = 2; % Always keep this consistent across mice
    data.target_stimchannel = 1; % Always keep this consistent across mice
    data.Fs = 10000; % Acquisition Frequency
    % Baseline
    data.baseline_winSize = [0 2];
    data.baseline_idx = [data.Fs*data.baseline_winSize(1)+1 data.Fs*data.baseline_winSize(2)];
    % Moving avg/smoothing
    data.smoothing_window = 1; % Size in milliseconds
    data.smoothing_window = data.smoothing_window * (data.Fs/1000);
    data.smoothing_overlay = 0.5; % Size in milliseconds
    data.smoothing_overlay = data.smoothing_overlay * (data.Fs/1000);
    % Peak/Dip detection thresholds
    data.max_point_within = [11 500];         % Max time within which peaks should occur (in ms)
    data.search_range = 2000;                 % Range to search peaks/dips within (in ms)
    % Depending on type of experiment, different things to process...
    if (data.mouse=="fig1c")
        data.onset_idx = nan(64,1); % Preallocation of the onset indicies
        data.offset_idx = nan(64,1); % Preallocation of the offset indicies
        for y = 1:5
            % Finds the number of stimuli, onset indicies, and offset indicies and saves it
            num_stim = 0;
            counter = 0;
            onset_idx = [];
            offset_idx = [];
            refSig_vector = data.acq(15*y-14).refSig{data.target_stimchannel,1};
                for z = 1:length(data.acq(15*y-14).refSig{data.target_stimchannel,1})
                    if (refSig_vector(z)>0.1 && counter == 0 || counter == 1)
                        if (counter == 0)
                            onset_idx = [onset_idx; z];
                        end
                        counter = 1;
                        if (refSig_vector(z)<0.1 && counter == 1)
                            num_stim = num_stim + 1;
                            offset_idx = [offset_idx; z-1]; % Offset index is one index before code detects return to baseline refSig levels
                            counter = 0;
                        end
                    end
                end
                data.num_stim(:,y) = [num_stim];
                data.onset_idx(:,y) = [onset_idx; nan(64-num_stim,1)];
                data.offset_idx(:,y) = [offset_idx; nan(64-num_stim,1)];
           % Finds the width of (each) pulse and saves it
           data.pulse_width(:,y) = nnz(data.acq(15*y-14).refSig{data.target_stimchannel,1}>0.1)/data.Fs;
        end
    elseif (data.mouse=="fig1d")
        for y = 1:4
            % Finds the number of stimuli, onset indicies, and offset indicies and saves it
            num_stim = 0;
            counter = 0;
            onset_idx = [];
            offset_idx = [];
            refSig_vector = data.acq(15*y-14).refSig{data.target_stimchannel,1};
                for z = 1:length(data.acq(15*y-14).refSig{data.target_stimchannel,1})
                    if (refSig_vector(z)>0.1 && counter == 0 || counter == 1)
                        if (counter == 0)
                            onset_idx = [onset_idx; z];
                        end
                        counter = 1;
                        if (refSig_vector(z)<0.1 && counter == 1)
                            num_stim = num_stim + 1;
                            offset_idx = [offset_idx; z-1]; % Offset index is one index before code detects return to baseline refSig levels
                            counter = 0;
                        end
                    end
                end
                data.num_stim(:,y) = num_stim;
                data.onset_idx(:,y) = onset_idx;
                data.offset_idx(:,y) = offset_idx;
           % Finds the width of (each) pulse and saves it
           data.pulse_width(:,y) = nnz(data.acq(15*y-14).refSig{data.target_stimchannel,1}>0.1)/data.Fs;
        end
    elseif (data.mouse=="fig1d-longpulse")
        for y = 1:3
            % Finds the number of stimuli, onset indicies, and offset indicies and saves it
            num_stim = 0;
            counter = 0;
            onset_idx = [];
            offset_idx = [];
            refSig_vector = data.acq(10*y-9).refSig{data.target_stimchannel,1};
                for z = 1:length(data.acq(10*y-9).refSig{data.target_stimchannel,1})
                    if (refSig_vector(z)>0.1 && counter == 0 || counter == 1)
                        if (counter == 0)
                            onset_idx = [onset_idx; z];
                        end
                        counter = 1;
                        if (refSig_vector(z)<0.1 && counter == 1)
                            num_stim = num_stim + 1;
                            offset_idx = [offset_idx; z-1]; % Offset index is one index before code detects return to baseline refSig levels
                            counter = 0;
                        end
                    end
                end
                data.num_stim(:,y) = num_stim;
                data.onset_idx(:,y) = onset_idx;
                data.offset_idx(:,y) = offset_idx;
           % Finds the width of (each) pulse and saves it
           data.pulse_width(:,y) = nnz(data.acq(10*y-9).refSig{data.target_stimchannel,1}>0.1)/data.Fs;
        end
    elseif (data.mouse=="fig1e")
        for y = 1:4
            % Finds the number of stimuli, onset indicies, and offset indicies and saves it
            num_stim = 0;
            counter = 0;
            onset_idx = [];
            offset_idx = [];
            refSig_vector = data.acq(15*y-14).refSig{data.target_stimchannel,1};
                for z = 1:length(data.acq(15*y-14).refSig{data.target_stimchannel,1})
                    if (refSig_vector(z)>0.1 && counter == 0 || counter == 1)
                        if (counter == 0)
                            onset_idx = [onset_idx; z];
                        end
                        counter = 1;
                        if (refSig_vector(z)<0.1 && counter == 1)
                            num_stim = num_stim + 1;
                            offset_idx = [offset_idx; z-1]; % Offset index is one index before code detects return to baseline refSig levels
                            counter = 0;
                        end
                    end
                end
                data.num_stim(:,y) = [num_stim];
                data.onset_idx(:,y) = [onset_idx; nan(4-num_stim,1)];
                data.offset_idx(:,y) = [offset_idx; nan(4-num_stim,1)];
           % Finds the width of (each) pulse and saves it
           data.pulse_width(:,y) = nnz(data.acq(15*y-14).refSig{data.target_stimchannel,1}>0.1)/data.Fs;
           % Pulse frequency (Onset or offset indicies and divide by the sampling freq to obtain the stim freq)
           if (data.num_stim==1)
               data.pulse_freq(:,y) = NaN;
           else
               data.pulse_freq(:,y) = (data.offset_idx(2,y)-data.offset_idx(1,y))/(data.Fs/1000);
           end
        end
    elseif (data.mouse=="fig1f")
        for y = 1:3
            % Finds the number of stimuli, onset indicies, and offset indicies and saves it
            num_stim = 0;
            counter = 0;
            onset_idx = [];
            offset_idx = [];
            refSig_vector = data.acq(15*y-14).refSig{data.target_stimchannel,1};
            for z = 1:length(data.acq(1).refSig{data.target_stimchannel,1})
                if (refSig_vector(z)>0.1 && counter == 0 || counter == 1)
                    if (counter == 0)
                        onset_idx = [onset_idx; z];
                    end
                    counter = 1;
                    if (refSig_vector(z)<0.1 && counter == 1)
                        num_stim = num_stim + 1;
                        offset_idx = [offset_idx; z-1]; % Offset index is one index before code detects return to baseline refSig levels
                        counter = 0;
                    end
                end
            end
            data.num_stim(:,y) = num_stim;
            data.onset_idx(:,y) = onset_idx;
            data.offset_idx(:,y) = offset_idx;
            % Paired-pulse interval (Use the onset or offset indicies)
            data.paired_pulse(:,y) = (data.onset_idx(2,y)-data.offset_idx(1,y)-1)/data.Fs;
        end
    end

    % For each sweep (corresponds to row in "data" struct), concatenate
    % the following data across all sweeps:
    % 1) FP signal (raw). DOES NOT clip. No preprocessing applied.
    % 2) FP ref signal (raw). Reference signal corresponding to stimuli applied
    if (isfield(data,'final'))
        data.final = [];
    end
    data.final.FP = [];
    data.final.refSig = [];
    data.final.allFP = [];
    data.final.allrefSig = [];
    if (data.mouse == "fig1c")
       process_idx = 5;
    elseif (data.mouse == "fig1d")
       process_idx = 4;
    elseif (data.mouse == "fig1d-longpulse")
       process_idx = 3;
    elseif (data.mouse == "fig1e")
       process_idx = 4;
    elseif (data.mouse == "fig1f")
       process_idx = 3;
    end

    for y = 1:length(data.acq)
         data.final(y).FP = data.acq(y).FP{data.target_FPchannel,1};
         data.final(y).refSig = data.acq(y).refSig{data.target_stimchannel,1};
    end

    for y = 1:process_idx
        if (data.mouse == "fig1d-longpulse")
            idx = 1:10;
            idx = idx+10*y-10;
            % Concatenation into 1 matrix
            data.final(y).allFP = [data.final(y).allFP data.final(idx).FP];
            data.final(y).allrefSig = [data.final(y).allrefSig data.final(idx).refSig];
        else
            idx = 1:15;
            idx = idx+15*y-15;
            % Concatenation into 1 matrix
            data.final(y).allFP = [data.final(y).allFP data.final(idx).FP];
            data.final(y).allrefSig = [data.final(y).allrefSig data.final(idx).refSig];
        end
    end

    for y = 1:process_idx
         % Calculate: 1) Mean, 2) Baselined FP, 3) %dF/F, and 4) SEM across all sweeps/trials
         % Mean Voltage (nanmean not required as raw data shouldn't have it)
         data.final(y).meanFP = mean(data.final(y).allFP,2);
         % Baselined FP Voltage (Calculate this after obtaining the meanFP)
         data.final(y).baselineFP = data.final(y).meanFP - mean(data.final(y).meanFP(data.baseline_idx(1):data.baseline_idx(2)));
         % SEM Voltage (Stays constant regardless of baselining or not)
         data.final(y).SEM_FP = SEM(data.final(y).allFP,2);
         % dF/F for each sweep. (dF/F = 100 x F(t) - F0/F0)
         for z = 1:size(data.final(1).allFP,2)
            data.final(y).allFP_dFF(:,z) = (data.final(y).allFP(:,z) - mean(data.final(y).allFP(data.baseline_idx(1):data.baseline_idx(2),z)))/mean(data.final(y).allFP(data.baseline_idx(1):data.baseline_idx(2),z));
         end
         % Mean dF/F (nanmean not required as raw data shouldn't have it)
         data.final(y).meanFP_dFF = mean(data.final(y).allFP_dFF,2);
         data.final(y).meanFP_dFF = data.final(y).meanFP_dFF*100;
         % Baselined dF/F
         data.final(y).baselineFP_dFF = data.final(y).meanFP_dFF - mean(data.final(y).meanFP_dFF(data.baseline_idx(1):data.baseline_idx(2)));
         % SEM dF/F
         data.final(y).SEM_FP_dFF = SEM(data.final(y).allFP_dFF,2)*100;
    end

    % PART 2: Calculating artifact properties
        if (data.mouse == "fig1c")
            for y = 1:process_idx
                for z = 1:data.num_stim(y)
                    if (data.num_stim(y) == z)
                        % PEAK
                        % AVERAGING/SMOOTHING
                        % Default: 1ms window size. 0.5ms overlay.
                        smoothed_interval = [];
                        % For the last stim, use the method of manually determining the interval of which to look at. (Default 2000ms interval)
                        extracted_interval = data.final(y).baselineFP_dFF(data.offset_idx(z,y):data.offset_idx(z,y)+data.search_range*(data.Fs/1000)-1);
                        for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(y).baselineFP_dFF(data.offset_idx(z,y):data.offset_idx(z,y)+data.search_range*(data.Fs/1000)-1))-data.smoothing_window
                            smoothed_interval(zz,1) = mean(extracted_interval(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                        end
                        % Align the Smoothed Interval to Real Intervals (i.e. the time from start of interval to time of 1st sampling) for precise indexing
                        data.final(y).smoothed_interval(:,z) = [smoothed_interval; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval),1)];
                        smoothed_interval_work = [zeros(round(data.smoothing_window/2),1); smoothed_interval(smoothed_interval~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval(smoothed_interval~=0))-data.smoothing_window/2,1)];
                        data.final(y).smoothed_interval_idx(:,z) = [find(smoothed_interval~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval(smoothed_interval~=0)),1)];
                        % Remove all values that are the stimulus itself (Approx Stimuli Interval * 3)
                        smoothed_interval_work(1:ceil(data.max_point_within(1)*(data.Fs/1000)/(data.smoothing_window-data.smoothing_overlay))) = 0;
                        % Remove all zeros (i.e. non-sampled indicies from the vector)
                        smoothed_interval_work(smoothed_interval_work == 0) = NaN;
                        % For the last stim, use the same method of manually determining the interval of which to look at.
                        % Counter prevents detecting smoothed artifacts as peaks
                        counter = 0; tolerance = 1;
                        while (counter == 0)
                            [data.final(y).smoothed_peak(:,z),data.final(y).smoothed_peak_idx(:,z)] = nanmax(smoothed_interval_work);
                            if (data.final(y).smoothed_peak_idx(:,z) > 1)
                                if abs(smoothed_interval_work(data.final(y).smoothed_peak_idx(:,z))-smoothed_interval_work(data.final(y).smoothed_peak_idx(:,z)-1)) > tolerance | find(find(data.final(y).smoothed_interval(:,z)==data.final(y).smoothed_peak(:,z))==data.final(y).smoothed_interval_idx(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2 > data.max_point_within(2)*(data.Fs/1000)
                                    smoothed_interval_work(data.final(y).smoothed_peak_idx(:,z)) = NaN;
                                else
                                    counter = 1;
                                end
                            else
                                counter = 1;
                            end
                        end 
                        counter = 0;
                        data.final(y).smoothed_peak_locs(:,z) = data.offset_idx(z,y) + find(find(data.final(y).smoothed_interval(:,z)==data.final(y).smoothed_peak(:,z))==data.final(y).smoothed_interval_idx(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;
                    
                        % DIP
                        % AVERAGING/SMOOTHING
                        smoothed_interval_post = [];
                        extracted_interval_post = data.final(y).baselineFP_dFF(data.offset_idx(z,y):data.offset_idx(z,y)+data.search_range*(data.Fs/1000)-1);
                        for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(y).baselineFP_dFF(data.offset_idx(z,y):data.offset_idx(z,y)+data.search_range*(data.Fs/1000)-1))-data.smoothing_window
                            smoothed_interval_post(zz,1) = mean(extracted_interval_post(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                        end
                        data.final(y).smoothed_interval_post(:,z) = [smoothed_interval_post; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_post),1)];
                        smoothed_interval_work_post = [zeros(round(data.smoothing_window/2),1); smoothed_interval_post(smoothed_interval_post~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_post(smoothed_interval_post~=0))-data.smoothing_window/2,1)];
                        data.final(y).smoothed_interval_idx_post(:,z) = [find(smoothed_interval_post~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_post(smoothed_interval_post~=0)),1)];
                        smoothed_interval_work_post(1:ceil(data.max_point_within(1)*(data.Fs/1000)/(data.smoothing_window-data.smoothing_overlay))) = 0;
                        smoothed_interval_work_post(smoothed_interval_work_post == 0) = NaN;
                        counter = 0; tolerance = 1;
                        while (counter == 0)
                            [data.final(y).smoothed_dip(:,z),data.final(y).smoothed_dip_idx(:,z)] = nanmin(smoothed_interval_work_post);
                            if abs(smoothed_interval_work_post(data.final(y).smoothed_dip_idx(:,z))-smoothed_interval_work_post(data.final(y).smoothed_dip_idx(:,z)-1)) > tolerance  | find(find(data.final(y).smoothed_interval_post(:,z)==data.final(y).smoothed_dip(:,z))==data.final(y).smoothed_interval_idx_post(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2 > data.max_point_within(2)*(data.Fs/1000)
                                smoothed_interval_work_post(data.final(y).smoothed_dip_idx(:,z)) = NaN;
                            else
                                counter = 1;
                            end
                        end 
                        counter = 0;
                        data.final(y).smoothed_dip_locs(:,z) = data.offset_idx(z,y) + find(find(data.final(y).smoothed_interval_post(:,z)==data.final(y).smoothed_dip(:,z))==data.final(y).smoothed_interval_idx_post(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;
                       
                        % DIP
                        % Smoothing an interval prior to the stimulation to acquire the baseline. (i.e. in comparision to the effect size of the pulse on the dip of the baseline fluorescence)
                        smoothed_interval_pre = [];
                        extracted_interval_pre = data.final(y).baselineFP_dFF(data.onset_idx(z,y)-20*(data.Fs/1000):data.onset_idx(z,y)-1);
                        for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(y).baselineFP_dFF(data.onset_idx(z,y)-20*(data.Fs/1000):data.onset_idx(z,y)-1))-data.smoothing_window
                            smoothed_interval_pre(zz,1) = mean(extracted_interval_pre(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                        end
                        data.final(y).smoothed_interval_pre(:,z) = [smoothed_interval_pre; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_pre),1)];
                        smoothed_interval_work_pre = [zeros(round(data.smoothing_window/2),1); smoothed_interval_pre(smoothed_interval_pre~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_pre(smoothed_interval_pre~=0))-data.smoothing_window/2,1)];
                        data.final(y).smoothed_interval_idx_pre(:,z) = [find(smoothed_interval_pre~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_pre(smoothed_interval_pre~=0)),1)];
                        smoothed_interval_work_pre(smoothed_interval_work_pre>20) = 0;
                        smoothed_interval_work_pre(smoothed_interval_work_pre<-20) = 0;
                        smoothed_interval_work_pre(smoothed_interval_work_pre == 0) = NaN;
                        data.final(y).smoothed_pre_stim(:,z) = nanmean(smoothed_interval_work_pre);
                        if (data.final(y).smoothed_pre_stim(:,z) - data.final(y).smoothed_dip(:,z) > 0)
                            data.final(y).dip_amplitude(:,z) = data.final(y).smoothed_pre_stim(:,z) - data.final(y).smoothed_dip(:,z);
                        else
                            data.final(y).dip_amplitude(:,z) = NaN; % If dip doesn't occur, than set its amplitude to 0.
                        end
                    else 
                       % PEAK
                       % AVERAGING/SMOOTHING
                       % Default: 1ms window size. 0.5ms overlay.
                       smoothed_interval = [];
                       % For the last stim, use the method of manually determining the interval of which to look at. (Default 100ms interval)
                       extracted_interval = data.final(y).baselineFP_dFF(data.offset_idx(z,y):data.onset_idx(z+1,y));
                       for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(1).baselineFP_dFF(data.offset_idx(z,y):data.onset_idx(z+1,y)))-data.smoothing_window
                           smoothed_interval(zz,1) = mean(extracted_interval(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                       end
                        % Align the Smoothed Interval to Real Intervals (i.e. the time from start of interval to time of 1st sampling) for precise indexing
                        data.final(y).smoothed_interval(:,z) = [smoothed_interval; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval),1)];
                        smoothed_interval_work = [zeros(round(data.smoothing_window/2),1); smoothed_interval(smoothed_interval~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval(smoothed_interval~=0))-data.smoothing_window/2,1)];
                        data.final(y).smoothed_interval_idx(:,z) = [find(smoothed_interval~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval(smoothed_interval~=0)),1)];
                        % Remove all values that are the stimulus itself (Approx Stimuli Interval * 3)
                        smoothed_interval_work(1:ceil(data.max_point_within(1)*(data.Fs/1000)/(data.smoothing_window-data.smoothing_overlay))) = 0;
                        % Remove all zeros (i.e. non-sampled indicies from the vector)
                        smoothed_interval_work(smoothed_interval_work == 0) = NaN;

                        % For the last stim, use the same method of manually determining the interval of which to look at.
                        % Counter prevents detecting smoothed artifacts as peaks
                        counter = 0; tolerance = 1;
                        while (counter == 0)
                            [data.final(y).smoothed_peak(:,z),data.final(y).smoothed_peak_idx(:,z)] = nanmax(smoothed_interval_work);
                            if (data.final(y).smoothed_peak_idx(:,z) > 1)
                                if abs(smoothed_interval_work(data.final(y).smoothed_peak_idx(:,z))-smoothed_interval_work(data.final(y).smoothed_peak_idx(:,z)-1)) > tolerance | find(find(data.final(y).smoothed_interval(:,z)==data.final(y).smoothed_peak(:,z))==data.final(y).smoothed_interval_idx(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2 > data.max_point_within(2)*(data.Fs/1000)
                                    smoothed_interval_work(data.final(y).smoothed_peak_idx(:,z)) = NaN;
                                else
                                    counter = 1;
                                end
                            else
                                counter = 1;
                            end
                        end
                        counter = 0;
                        data.final(y).smoothed_peak_locs(:,z) = data.offset_idx(z,y) + find(find(data.final(y).smoothed_interval(:,z)==data.final(y).smoothed_peak(:,z))==data.final(y).smoothed_interval_idx(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;

                       % DIP
                       % AVERAGING/SMOOTHING
                       smoothed_interval_post = [];
                       extracted_interval_post = data.final(y).baselineFP_dFF(data.offset_idx(z,y):data.onset_idx(z+1,y));
                       for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(1).baselineFP_dFF(data.offset_idx(z,y):data.onset_idx(z+1,y)))-data.smoothing_window
                           smoothed_interval_post(zz,1) = mean(extracted_interval_post(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                       end
                       data.final(y).smoothed_interval_post(:,z) = [smoothed_interval_post; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_post),1)];
                       smoothed_interval_work_post = [zeros(round(data.smoothing_window/2),1); smoothed_interval_post(smoothed_interval_post~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_post(smoothed_interval_post~=0))-data.smoothing_window/2,1)];
                       data.final(y).smoothed_interval_idx_post(:,z) = [find(smoothed_interval_post~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_post(smoothed_interval_post~=0)),1)];
                       smoothed_interval_work_post(1:ceil(data.max_point_within(1)*(data.Fs/1000)/(data.smoothing_window-data.smoothing_overlay))) = 0;
                       smoothed_interval_work_post(smoothed_interval_work_post == 0) = NaN;
                       counter = 0; tolerance = 1;
                       while (counter == 0)
                           [data.final(y).smoothed_dip(:,z),data.final(y).smoothed_dip_idx(:,z)] = nanmin(smoothed_interval_work_post);
                           if abs(smoothed_interval_work_post(data.final(y).smoothed_dip_idx(:,z))-smoothed_interval_work_post(data.final(y).smoothed_dip_idx(:,z)-1)) > tolerance | find(find(data.final(y).smoothed_interval_post(:,z)==data.final(y).smoothed_dip(:,z))==data.final(y).smoothed_interval_idx_post(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2 > data.max_point_within(2)*(data.Fs/1000)
                               smoothed_interval_work_post(data.final(y).smoothed_dip_idx(:,z)) = NaN;
                           else
                               counter = 1;
                           end
                       end 
                       counter = 0;
                       data.final(y).smoothed_dip_locs(:,z) = data.offset_idx(z,y) + find(find(data.final(y).smoothed_interval_post(:,z)==data.final(y).smoothed_dip(:,z))==data.final(y).smoothed_interval_idx_post(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;
                       
                       % DIP
                       % Smoothing an interval prior to the stimulation to acquire the baseline. (i.e. in comparision to the effect size of the pulse on the dip of the baseline fluorescence)
                       smoothed_interval_pre = [];
                       extracted_interval_pre = data.final(y).baselineFP_dFF(data.onset_idx(z,y)-20*(data.Fs/1000):data.onset_idx(z,y)-1);
                       for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(y).baselineFP_dFF(data.onset_idx(z,y)-20*(data.Fs/1000):data.onset_idx(z,y)-1))-data.smoothing_window
                           smoothed_interval_pre(zz,1) = mean(extracted_interval_pre(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                       end
                       data.final(y).smoothed_interval_pre(:,z) = [smoothed_interval_pre; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_pre),1)];
                       smoothed_interval_work_pre = [zeros(round(data.smoothing_window/2),1); smoothed_interval_pre(smoothed_interval_pre~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_pre(smoothed_interval_pre~=0))-data.smoothing_window/2,1)];
                       data.final(y).smoothed_interval_idx_pre(:,z) = [find(smoothed_interval_pre~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_pre(smoothed_interval_pre~=0)),1)];
                       smoothed_interval_work_pre(smoothed_interval_work_pre>20) = 0;
                       smoothed_interval_work_pre(smoothed_interval_work_pre<-20) = 0;
                      smoothed_interval_work_pre(smoothed_interval_work_pre == 0) = NaN;
                       data.final(y).smoothed_pre_stim(:,z) = nanmean(smoothed_interval_work_pre);
                       if (data.final(y).smoothed_pre_stim(:,z) - data.final(y).smoothed_dip(:,z) > 0)
                           data.final(y).dip_amplitude(:,z) = data.final(y).smoothed_pre_stim(:,z) - data.final(y).smoothed_dip(:,z);
                       else
                           data.final(y).dip_amplitude(:,z) = NaN; % If dip doesn't occur, than set its amplitude to NaN.
                       end
                    end
                end

                % Calculate the 10% and 90% of the peak amplitude.
                % Calculate the Decay Time Constant (37 of peak amplitude).
                data.final(y).max_peak_fraction = [0.1*max(data.final(y).smoothed_peak) 0.37*max(data.final(y).smoothed_peak) 0.9*max(data.final(y).smoothed_peak)];
                % 10%->90% Signal Rise Time: This is the time corresponding to how fast a rise in signal.
                % Get the indicies between the offset of the stimulation to the est. peak amplitude.
                % Intersection at 10% dF/F of peak amp when the signal is rising towards peak.
                counter1 = 0; idx_counter= 1; index = 0;
                nan_counter = 0; number_stim = 0; % keeps track of how many stimuli we have gone through
                % Tolerance value: 0.020 (%dF/F). If the value falls within this range
                tolerance = 0.015;
                % Prevents detection of the stimuli as peaks (which results in negative rise/decay time values)
                adjacent_noise_tolerance = 2;
                while (counter1 == 0 && idx_counter < size(data.final(y).smoothed_interval,1)*size(data.final(y).smoothed_interval,2))
                    if (data.final(y).smoothed_interval(idx_counter) > 0 && data.final(y).smoothed_interval(idx_counter)~=0 && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).max_peak_fraction(1,1))<=tolerance && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).smoothed_interval(idx_counter+data.smoothing_window/2))<adjacent_noise_tolerance)
                        rise_to_peak = index;
                        counter1 = 1; % Break the loop when signal-to-0 decay idx has been identified
                    elseif (data.final(y).smoothed_interval(idx_counter) > data.final(y).max_peak_fraction(1,1) && data.final(y).smoothed_interval(idx_counter)~=0 && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).max_peak_fraction(1,1))>tolerance && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).smoothed_interval(idx_counter+data.smoothing_window/2))<adjacent_noise_tolerance)
                        rise_to_peak = index;
                        counter1 = 1; % Break the loop when signal-to-0 decay idx has been identified
                    else
                        if (data.final(y).smoothed_interval(idx_counter)~=0 && ~isnan(data.final(y).smoothed_interval(idx_counter)))
                            index = index + 1;
                        elseif (isnan(data.final(y).smoothed_interval(idx_counter)) && nan_counter == 0)
                            number_stim = number_stim + 1;
                            nan_counter = 1;
                        elseif (~isnan(data.final(y).smoothed_interval(idx_counter)) && nan_counter == 1)
                            nan_counter = 0;
                        end
                        idx_counter = idx_counter + 1;
                    end
                end
                if (idx_counter > size(data.final(y).smoothed_interval,1)*size(data.final(y).smoothed_interval,2) || data.final(y).max_peak_fraction(1,1) < 0.1)
                   rise_intersection_at_10_idx = NaN;
                else
                   rise_intersection_at_10_idx = data.offset_idx(1,y) + rise_to_peak * (data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;
                end

                % Intersection at 90% dF/F of peak amp when the signal is rising towards peak.
                counter1 = 0; idx_counter= 1; index = 0;
                nan_counter = 0; number_stim = 0;                       % keeps track of how many stimuli we have gone through
                tolerance = 0.015;                                      % Tolerance value: (%dF/F). If the value falls within this range
                adjacent_noise_tolerance = 2;                            % Prevents detection of the stimuli as peaks (which results in negative rise/decay time values)
                while (counter1 == 0 && idx_counter < size(data.final(y).smoothed_interval,1)*size(data.final(y).smoothed_interval,2))
                    if (data.final(y).smoothed_interval(idx_counter) > 0 && data.final(y).smoothed_interval(idx_counter)~=0 && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).max_peak_fraction(1,3))<=tolerance &&  abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).smoothed_interval(idx_counter+data.smoothing_window/2))<adjacent_noise_tolerance)
                        rise_to_peak = index;
                        counter1 = 1; % Break the loop when signal-to-0 decay idx has been identified
                    elseif (data.final(y).smoothed_interval(idx_counter) > data.final(y).max_peak_fraction(1,3) && data.final(y).smoothed_interval(idx_counter)~=0 && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).max_peak_fraction(1,3))>tolerance && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).smoothed_interval(idx_counter+data.smoothing_window/2))<adjacent_noise_tolerance)
                        rise_to_peak = index;
                        counter1 = 1; % Break the loop when signal-to-0 decay idx has been identified
                    else
                        if (data.final(y).smoothed_interval(idx_counter)~=0 && ~isnan(data.final(y).smoothed_interval(idx_counter)))
                            index = index + 1;
                        elseif (isnan(data.final(y).smoothed_interval(idx_counter)) && nan_counter == 0)
                            number_stim = number_stim + 1;
                            nan_counter = 1;
                        elseif (~isnan(data.final(y).smoothed_interval(idx_counter)) && nan_counter == 1)
                            nan_counter = 0;
                        end
                        idx_counter = idx_counter + 1;
                    end
                end
                % Calculate 10%->90% Signal rise time
                if (idx_counter > size(data.final(y).smoothed_interval,1)*size(data.final(y).smoothed_interval,2) || data.final(y).max_peak_fraction(1,3) < 0.9)
                   rise_intersection_at_90_idx = NaN;
                else
                rise_intersection_at_90_idx = data.offset_idx(1,y) + number_stim * data.pulse_width(y) * data.Fs +  rise_to_peak * (data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;
                end
                data.final(y).ten_to_ninety_rise_time = (rise_intersection_at_90_idx - rise_intersection_at_10_idx)/data.Fs;
                
                % Decay Time Constant: (This is approximately when the signal reaches 37% from its max value)
                % Calculate the fall intersection at which we see signal baselining again after peak (I.e. Intersection at 0% dF/F)
                counter2 = 0;                          
                % idx_counter = 0;                                    % Don't reset the conuter, as it will be continued for finding the 37% value post-peak
                tolerance = 0.015;                                    % Tolerance value:= (%dF/F). If the value falls within this range near 0, than signal has reached baseline, post peak.
                idx_counter = find(data.final(y).smoothed_peak(end)==data.final(y).smoothed_interval); % Always start index counter from the last peak that was seen.
                adjacent_noise_tolerance = 2;                     % Prevents detection of the stimuli as peaks (which results in negative rise/decay time values)
                while (counter2 == 0 && idx_counter < size(data.final(y).smoothed_interval,1)*size(data.final(y).smoothed_interval,2))
                    if (data.final(y).smoothed_interval(idx_counter) > 0 && data.final(y).smoothed_interval(idx_counter)~=0 && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).max_peak_fraction(1,2))<=tolerance && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).smoothed_interval(idx_counter+data.smoothing_window/2))<adjacent_noise_tolerance)
                        decay_to_37 = index;
                        counter2 = 1; % Break the loop when signal-to-0 decay idx has been identified
                    elseif (data.final(y).smoothed_interval(idx_counter) < data.final(y).max_peak_fraction(1,2) && data.final(y).smoothed_interval(idx_counter)~=0 && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).max_peak_fraction(1,2))>tolerance && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).smoothed_interval(idx_counter+data.smoothing_window/2))<adjacent_noise_tolerance)    
                        decay_to_37 = index;
                        counter2 = 1; % Break the loop when signal-to-0 decay idx has been identified
                    else
                        if (data.final(y).smoothed_interval(idx_counter)~=0 && ~isnan(data.final(y).smoothed_interval(idx_counter)))
                            index = index + 1;
                        elseif (isnan(data.final(y).smoothed_interval(idx_counter)) && nan_counter == 0)
                            number_stim = number_stim + 1;
                            nan_counter = 1;
                        elseif (~isnan(data.final(y).smoothed_interval(idx_counter)) && nan_counter == 1)
                            nan_counter = 0;
                        end
                        idx_counter = idx_counter + 1;
                    end
                end
                % PEAK
                if (idx_counter < size(data.final(y).smoothed_interval,1)*size(data.final(y).smoothed_interval,2))
                    decay_intersection_at_37_idx = data.offset_idx(find(~isnan(data.offset_idx(:,y)),1,'last'),y) + number_stim * data.pulse_width(y) * data.Fs + decay_to_37 * (data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;
                    % Calculate the Decay Time Constant (at 37% of peak signal)
                    if (~isempty(decay_to_37) & (decay_intersection_at_37_idx - data.final(y).smoothed_peak_locs(end))/data.Fs > 0 & data.final(y).max_peak_fraction(1,2) > 0.37)
                         data.final(y).decay_time_constant = (decay_intersection_at_37_idx - data.final(y).smoothed_peak_locs(end))/data.Fs;
                    else
                        data.final(y).decay_time_constant = NaN;
                    end
                else
                    data.final(y).decay_time_constant = NaN;
                end
                counter2 = 0; idx_counter = 0; index = 0;   % Reset counters
                % DIP
                % Calculate the max dip amplitude (relative to each stim) across all
                data.final(y).max_dip_amplitude = nanmax(data.final(y).dip_amplitude);
                if data.final(y).max_dip_amplitude < 0
                    data.final(y).max_dip_amplitude = NaN; % NaN if no dips occur
                end
                % Calculate the mean dip amplitude per each condition
                data.final(y).mean_dip_amplitude = nanmean(data.final(y).dip_amplitude);
                % Calculate the total (i.e. accumulated dip) across all 
                data.final(y).total_dip = data.final(y).smoothed_pre_stim(1)-min(data.final(y).smoothed_dip);
                if data.final(y).total_dip < 0
                    data.final(y).total_dip = NaN; % NaN if no dips occur
                end
            end
        elseif (data.mouse == "fig1d")
          for y = 1:process_idx
                % PEAK
                % AVERAGING/SMOOTHING
                % Default: 1ms window size. 0.5ms overlay.
                smoothed_interval = [];
                % For the last stim, use the method of manually determining the interval of which to look at. (Default 2000ms interval)
                extracted_interval = data.final(y).baselineFP_dFF(data.offset_idx(1,y):data.offset_idx(1,y)+data.search_range*(data.Fs/1000)-1);
                for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(y).baselineFP_dFF(data.offset_idx(1,y):data.offset_idx(1,y)+data.search_range*(data.Fs/1000)-1))-data.smoothing_window
                    smoothed_interval(zz,1) = mean(extracted_interval(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                end
                % Align the Smoothed Interval to Real Intervals (i.e. the time from start of interval to time of 1st sampling) for precise indexing
                data.final(y).smoothed_interval(:,1) = [smoothed_interval; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval),1)];
                smoothed_interval_work = [zeros(round(data.smoothing_window/2),1); smoothed_interval(smoothed_interval~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval(smoothed_interval~=0))-data.smoothing_window/2,1)];
                data.final(y).smoothed_interval_idx(:,1) = [find(smoothed_interval~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval(smoothed_interval~=0)),1)];
                % Remove all values that are the stimulus itself (Approx Stimuli Interval * 3)
                smoothed_interval_work(1:ceil(data.max_point_within(1)*(data.Fs/1000)/(data.smoothing_window-data.smoothing_overlay))) = 0;
                % Remove all zeros (i.e. non-sampled indicies from the vector)
                smoothed_interval_work(smoothed_interval_work == 0) = NaN;

                % For the last stim, use the same method of manually determining the interval of which to look at.
                % Counter prevents detecting smoothed artifacts as peaks
                counter = 0; tolerance = 1;
                while (counter == 0)
                    [data.final(y).smoothed_peak,data.final(y).smoothed_peak_idx] = nanmax(smoothed_interval_work);
                    if (data.final(y).smoothed_peak_idx > 1)
                        if abs(smoothed_interval_work(data.final(y).smoothed_peak_idx)-smoothed_interval_work(data.final(y).smoothed_peak_idx-1)) > tolerance | find(find(data.final(y).smoothed_interval==data.final(y).smoothed_peak)==data.final(y).smoothed_interval_idx)*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2 > data.max_point_within(2)*(data.Fs/1000)
                            smoothed_interval_work(data.final(y).smoothed_peak_idx) = NaN;
                        else
                            counter = 1;
                        end
                    else
                        counter = 1;
                    end
                end
                counter = 0;
                data.final(y).smoothed_peak_locs = data.offset_idx(1,y) + find(find(data.final(y).smoothed_interval==data.final(y).smoothed_peak)==data.final(y).smoothed_interval_idx)*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;

                % DIP
                % AVERAGING/SMOOTHING
                smoothed_interval_post = [];
                extracted_interval_post = data.final(y).baselineFP_dFF(data.offset_idx(1,y):data.offset_idx(1,y)+data.search_range*(data.Fs/1000)-1);
                for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(y).baselineFP_dFF(data.offset_idx(1,y):data.offset_idx(1,y)+data.search_range*(data.Fs/1000)-1))-data.smoothing_window
                    smoothed_interval_post(zz,1) = mean(extracted_interval_post(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                end
                data.final(y).smoothed_interval_post(:,1) = [smoothed_interval_post; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_post),1)];
                smoothed_interval_work_post = [zeros(round(data.smoothing_window/2),1); smoothed_interval_post(smoothed_interval_post~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_post(smoothed_interval_post~=0))-data.smoothing_window/2,1)];
                data.final(y).smoothed_interval_idx_post(:,1) = [find(smoothed_interval_post~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_post(smoothed_interval_post~=0)),1)];
                smoothed_interval_work_post(1:ceil(data.max_point_within(1)*(data.Fs/1000)/(data.smoothing_window-data.smoothing_overlay))) = 0;
                smoothed_interval_work_post(smoothed_interval_work_post == 0) = NaN;
                counter = 0; tolerance = 1;
                while (counter == 0)
                    [data.final(y).smoothed_dip,data.final(y).smoothed_dip_idx] = nanmin(smoothed_interval_work_post);
                    if abs(smoothed_interval_work_post(data.final(y).smoothed_dip_idx)-smoothed_interval_work_post(data.final(y).smoothed_dip_idx-1)) > tolerance  | find(find(data.final(y).smoothed_interval_post==data.final(y).smoothed_dip)==data.final(y).smoothed_interval_idx_post)*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2 > data.max_point_within(2)*(data.Fs/1000)
                        smoothed_interval_work_post(data.final(y).smoothed_dip_idx) = NaN;
                    else
                        counter = 1;
                    end
                end 
                counter = 0;
                data.final(y).smoothed_dip_locs(:,1) = data.offset_idx(1,y) + find(find(data.final(y).smoothed_interval_post(:,1)==data.final(y).smoothed_dip(:,1))==data.final(y).smoothed_interval_idx_post(:,1))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;
                
                % DIP
                % Smoothing an interval prior to the stimulation to acquire the baseline. (i.e. in comparision to the effect size of the pulse on the dip of the baseline fluorescence)
                smoothed_interval_pre = [];
                extracted_interval_pre = data.final(y).baselineFP_dFF(data.offset_idx(1,y)-20*(data.Fs/1000):data.offset_idx(1,y)-1);
                for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(y).baselineFP_dFF(data.offset_idx(1,y)-20*(data.Fs/1000):data.offset_idx(1,y)-1))-data.smoothing_window
                    smoothed_interval_pre(zz,1) = mean(extracted_interval_pre(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                end
                data.final(y).smoothed_interval_pre(:,1) = [smoothed_interval_pre; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_pre),1)];
                smoothed_interval_work_pre = [zeros(round(data.smoothing_window/2),1); smoothed_interval_pre(smoothed_interval_pre~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_pre(smoothed_interval_pre~=0))-data.smoothing_window/2,1)];
                data.final(y).smoothed_interval_idx_pre(:,1) = [find(smoothed_interval_pre~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_pre(smoothed_interval_pre~=0)),1)];
                smoothed_interval_work_pre(smoothed_interval_work_pre>20) = 0;
                smoothed_interval_work_pre(smoothed_interval_work_pre<-20) = 0;
                smoothed_interval_work_pre(smoothed_interval_work_pre == 0) = NaN;
                data.final(y).smoothed_pre_stim(:,z) = nanmean(smoothed_interval_work_pre);

                % Calculate the Decay Time Constant (37 of peak amplitude).
                data.final(y).max_peak_fraction = [0.1*max(data.final(y).smoothed_peak) 0.37*max(data.final(y).smoothed_peak) 0.9*max(data.final(y).smoothed_peak)];
                % Intersection at 10% dF/F of peak amp when the signal is rising towards peak.
                rise_to_peak = []; counter1 = 0; idx_counter= 1; index = 0;
                nan_counter = 0; number_stim = 0; % keeps track of how many stimuli we have gone through
                % Tolerance value: (%dF/F). If the value falls within this range
                tolerance = 0.015;
                % Prevents detection of the stimuli as peaks (which results in negative rise/decay time values)
                adjacent_noise_tolerance = 2;
                while (counter1 == 0 && idx_counter < size(data.final(y).smoothed_interval,1)*size(data.final(y).smoothed_interval,2))
                    if (data.final(y).smoothed_interval(idx_counter) > 0 && data.final(y).smoothed_interval(idx_counter)~=0 && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).max_peak_fraction(1,1))<=tolerance && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).smoothed_interval(idx_counter+data.smoothing_window/2))<adjacent_noise_tolerance)
                        rise_to_peak = index;
                        counter1 = 1; % Break the loop when signal-to-0 decay idx has been identified
                    elseif (data.final(y).smoothed_interval(idx_counter) > data.final(y).max_peak_fraction(1,1) && data.final(y).smoothed_interval(idx_counter)~=0 && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).max_peak_fraction(1,1))>tolerance && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).smoothed_interval(idx_counter+data.smoothing_window/2))<adjacent_noise_tolerance)
                        rise_to_peak = index;
                        counter1 = 1; % Break the loop when signal-to-0 decay idx has been identified
                    else
                        if (data.final(y).smoothed_interval(idx_counter)~=0 && ~isnan(data.final(y).smoothed_interval(idx_counter)))
                            index = index + 1;
                        elseif (isnan(data.final(y).smoothed_interval(idx_counter)) && nan_counter == 0)
                            number_stim = number_stim + 1;
                            nan_counter = 1;
                        elseif (~isnan(data.final(y).smoothed_interval(idx_counter)) && nan_counter == 1)
                            nan_counter = 0;
                        end
                        idx_counter = idx_counter + 1;
                    end
                end
                
                if (idx_counter > size(data.final(y).smoothed_interval,1)*size(data.final(y).smoothed_interval,2) || data.final(y).max_peak_fraction(1,1) < 0.1)
                   rise_intersection_at_10_idx = NaN;
                else
                   rise_intersection_at_10_idx = data.offset_idx(1,y) + rise_to_peak * (data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;
                end

                % Intersection at 90% dF/F of peak amp when the signal is rising towards peak.
                rise_to_peak = []; counter1 = 0; idx_counter= 1; index = 0;
                nan_counter = 0; number_stim = 0; % keeps track of how many stimuli we have gone through
                % Tolerance value: (%dF/F). If the value falls within this range
                tolerance = 0.015;
                adjacent_noise_tolerance = 2;
                while (counter1 == 0 && idx_counter < size(data.final(y).smoothed_interval,1)*size(data.final(y).smoothed_interval,2))
                    if (data.final(y).smoothed_interval(idx_counter) > 0 && data.final(y).smoothed_interval(idx_counter)~=0 && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).max_peak_fraction(1,3))<=tolerance && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).smoothed_interval(idx_counter+data.smoothing_window/2))<adjacent_noise_tolerance)
                        rise_to_peak = index;
                        counter1 = 1; % Break the loop when signal-to-0 decay idx has been identified
                    elseif (data.final(y).smoothed_interval(idx_counter) > data.final(y).max_peak_fraction(1,3) && data.final(y).smoothed_interval(idx_counter)~=0 && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).max_peak_fraction(1,3))>tolerance && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).smoothed_interval(idx_counter+data.smoothing_window/2))<adjacent_noise_tolerance)
                        rise_to_peak = index;
                        counter1 = 1; % Break the loop when signal-to-0 decay idx has been identified
                    else
                        if (data.final(y).smoothed_interval(idx_counter)~=0 && ~isnan(data.final(y).smoothed_interval(idx_counter)))
                            index = index + 1;
                        elseif (isnan(data.final(y).smoothed_interval(idx_counter)) && nan_counter == 0)
                            number_stim = number_stim + 1;
                            nan_counter = 1;
                        elseif (~isnan(data.final(y).smoothed_interval(idx_counter)) && nan_counter == 1)
                            nan_counter = 0;
                        end
                        idx_counter = idx_counter + 1;
                    end
                end

                if (idx_counter > size(data.final(y).smoothed_interval,1)*size(data.final(y).smoothed_interval,2) || data.final(y).max_peak_fraction(1,3) < 0.9)
                   rise_intersection_at_90_idx = NaN;
                else
                rise_intersection_at_90_idx = data.offset_idx(1,y) + number_stim * data.pulse_width(y) * data.Fs +  rise_to_peak * (data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;
                end
                % Calculate 10%->90% Signal rise time
                data.final(y).ten_to_ninety_rise_time = (rise_intersection_at_90_idx - rise_intersection_at_10_idx)/data.Fs;

                % Decay Time Constant: (This is approximately when the signal reaches 37% from its max value)
                % Calculate the fall intersection at which we see signal baselining again after peak
                % (I.e. Intersection at 0% dF/F)
                % Don't reset the conuter, as it will be continued for finding the 37% value post-peak
                decay_to_37 = []; counter2 = 0;
                tolerance = 0.015;                                                                          % Tolerance value: (%dF/F). If the value falls within this range 
                idx_counter = find(data.final(y).smoothed_peak(end)==data.final(y).smoothed_interval);      % Always start index counter from the last peak that was seen.
                adjacent_noise_tolerance = 2;                                                               % Prevents detection of the stimuli as peaks (which results in negative rise/decay time values)
                while (counter2 == 0 && idx_counter < size(data.final(y).smoothed_interval,1)*size(data.final(y).smoothed_interval,2))
                    if (data.final(y).smoothed_interval(idx_counter) > 0 && data.final(y).smoothed_interval(idx_counter)~=0 && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).max_peak_fraction(1,2))<=tolerance && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).smoothed_interval(idx_counter+data.smoothing_window/2))<adjacent_noise_tolerance)
                        decay_to_37 = index;
                        counter2 = 1; % Break the loop when signal-to-0 decay idx has been identified
                    elseif (data.final(y).smoothed_interval(idx_counter) < data.final(y).max_peak_fraction(1,2) && data.final(y).smoothed_interval(idx_counter)~=0 && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).max_peak_fraction(1,2))>tolerance && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).smoothed_interval(idx_counter+data.smoothing_window/2))<adjacent_noise_tolerance)    
                        decay_to_37 = index;
                        counter2 = 1; % Break the loop when signal-to-0 decay idx has been identified
                    else
                        if (data.final(y).smoothed_interval(idx_counter)~=0 && ~isnan(data.final(y).smoothed_interval(idx_counter)))
                            index = index + 1;
                        elseif (isnan(data.final(y).smoothed_interval(idx_counter)) && nan_counter == 0)
                            number_stim = number_stim + 1;
                            nan_counter = 1;
                        elseif (~isnan(data.final(y).smoothed_interval(idx_counter)) && nan_counter == 1)
                            nan_counter = 0;
                        end
                        idx_counter = idx_counter + 1;
                    end
                end

                if (idx_counter < size(data.final(y).smoothed_interval,1)*size(data.final(y).smoothed_interval,2))
                    decay_intersection_at_37_idx = data.offset_idx(find(~isnan(data.offset_idx(:,y)),1,'last'),y) + number_stim * data.pulse_width(y) * data.Fs + decay_to_37 * (data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;
                    % Calculate the Decay Time Constant (at 37% of peak signal)
                    if (~isempty(decay_to_37) & (decay_intersection_at_37_idx - data.final(y).smoothed_peak_locs(end))/data.Fs > 0 & data.final(y).max_peak_fraction(1,2) > 0.37)
                         data.final(y).decay_time_constant = (decay_intersection_at_37_idx - data.final(y).smoothed_peak_locs(end))/data.Fs;
                    else
                        data.final(y).decay_time_constant = NaN;
                    end
                else
                    data.final(y).decay_time_constant = NaN;
                end
                counter2 = 0; idx_counter = 0; index = 0;   % Reset counters

                % DIP
                if (data.final(y).smoothed_pre_stim(:,1) - data.final(y).smoothed_dip(:,1) > 0)
                    data.final(y).dip_amplitude(:,1) = data.final(y).smoothed_pre_stim(:,1) - data.final(y).smoothed_dip(:,1);
                else
                    data.final(y).dip_amplitude(:,1) = NaN; % If dip doesn't occur, than set its amplitude to 0.
                end

                % Calculate the max dip amplitude (relative to each stim) across all
                data.final(y).max_dip_amplitude = nanmax(data.final(y).dip_amplitude);
                if data.final(y).max_dip_amplitude < 0
                    data.final(y).max_dip_amplitude = NaN; % NaN if no dips occur
                end

                % Calculate the mean dip amplitude per each condition
                data.final(y).mean_dip_amplitude = nanmean(data.final(y).dip_amplitude);

                % Calculate the total (i.e. accumulated dip) across all 
                data.final(y).total_dip = data.final(y).smoothed_pre_stim(1)-min(data.final(y).smoothed_dip);
                if data.final(y).total_dip < 0
                    data.final(y).total_dip = NaN; % NaN if no dips occur
                end

          end

        elseif (data.mouse == "fig1d-longpulse")
          for y = 1:process_idx
                % AVERAGING/SMOOTHING
                % Default: 1ms window size. 0.5ms overlay.
                smoothed_interval = [];
                % For the last stim, use the method of manually determining the interval of which to look at. (Default 2000ms interval)
                extracted_interval = data.final(y).baselineFP_dFF(data.offset_idx(1,y):data.offset_idx(1,y)+data.search_range*(data.Fs/1000)-1);
                for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(y).baselineFP_dFF(data.offset_idx(1,y):data.offset_idx(1,y)+data.search_range*(data.Fs/1000)-1))-data.smoothing_window
                    smoothed_interval(zz,1) = mean(extracted_interval(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                end
                % Align the Smoothed Interval to Real Intervals (i.e. the time from start of interval to time of 1st sampling) for precise indexing
                data.final(y).smoothed_interval(:,1) = [smoothed_interval; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval),1)];
                smoothed_interval_work = [zeros(round(data.smoothing_window/2),1); smoothed_interval(smoothed_interval~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval(smoothed_interval~=0))-data.smoothing_window/2,1)];
                data.final(y).smoothed_interval_idx(:,1) = [find(smoothed_interval~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval(smoothed_interval~=0)),1)];
                % Remove all values that are the stimulus itself (Approx Stimuli Interval * 3)
                smoothed_interval_work(1:ceil(data.max_point_within(1)*(data.Fs/1000)/(data.smoothing_window-data.smoothing_overlay))) = 0;
                % Remove all values related with the stimulus fluctuations (i.e. very high numbers and low numbers)
                smoothed_interval_work(smoothed_interval_work>40) = 0;
                smoothed_interval_work(smoothed_interval_work<-5) = 0;
                % Remove all zeros (i.e. non-sampled indicies from the vector)
                smoothed_interval_work(smoothed_interval_work == 0) = NaN;
                % For the last stim, use the same method of manually determining the interval of which to look at.
                % Counter prevents detecting smoothed artifacts as peaks
                counter = 0; tolerance = 1;
                while (counter == 0)
                    [data.final(y).smoothed_peak,data.final(y).smoothed_peak_idx] = nanmax(smoothed_interval_work);
                    if (data.final(y).smoothed_peak_idx > 1)
                        if abs(smoothed_interval_work(data.final(y).smoothed_peak_idx)-smoothed_interval_work(data.final(y).smoothed_peak_idx-1)) > tolerance | find(find(data.final(y).smoothed_interval==data.final(y).smoothed_peak)==data.final(y).smoothed_interval_idx)*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2 > data.max_point_within(2)*(data.Fs/1000)
                            smoothed_interval_work(data.final(y).smoothed_peak_idx) = NaN;
                        else
                            counter = 1;
                        end
                    else
                        counter = 1;
                    end
                end 
                if (data.final(y).smoothed_peak_idx < (data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2)
                    % If no peak meaning peak location is smaller than smoothing interval (i.e. 10th index), use 4ms after the stim offset as the peak location and value
                    data.final(y).smoothed_peak_locs = data.offset_idx(1,y) + 40;
                    data.final(y).smoothed_peak = data.final(y).baselineFP_dFF(data.final(y).smoothed_peak_locs);
                else
                    data.final(y).smoothed_peak_locs = data.offset_idx(1,y) + find(find(data.final(y).smoothed_interval==data.final(y).smoothed_peak)==data.final(y).smoothed_interval_idx)*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;
                end

                data.final(y).max_peak_fraction = [0.1*max(data.final(y).smoothed_peak) 0.37*max(data.final(y).smoothed_peak) 0.9*max(data.final(y).smoothed_peak)];
                % 10%->90% Signal Rise Time: This is the time corresponding to how fast a rise in signal.
                % Get the indicies between the offset of the stimulation to the est. peak amplitude.

                % Intersection at 10% dF/F of peak amp when the signal is rising towards peak.
                rise_to_peak = []; counter1 = 0; idx_counter= 1; index = 0; 
                nan_counter = 0; number_stim = 0;                               % keeps track of how many stimuli we have gone through
                tolerance = 0.015;                                              % Tolerance value: (%dF/F). If the value falls within this range
                adjacent_noise_tolerance = 2;                                   % Prevents detection of the stimuli as peaks (which results in negative rise/decay time values)
                while (counter1 == 0 && idx_counter < size(data.final(y).smoothed_interval,1)*size(data.final(y).smoothed_interval,2))
                    if (data.final(y).smoothed_interval(idx_counter) > 0 && data.final(y).smoothed_interval(idx_counter)~=0 && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).max_peak_fraction(1,1))<=tolerance && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).smoothed_interval(idx_counter+data.smoothing_window/2))<adjacent_noise_tolerance)
                        rise_to_peak = index;
                        counter1 = 1; % Break the loop when signal-to-0 decay idx has been identified
                    elseif (data.final(y).smoothed_interval(idx_counter) > data.final(y).max_peak_fraction(1,1) && data.final(y).smoothed_interval(idx_counter)~=0 && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).max_peak_fraction(1,1))>tolerance && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).smoothed_interval(idx_counter+data.smoothing_window/2))<adjacent_noise_tolerance)
                        rise_to_peak = index;
                        counter1 = 1; % Break the loop when signal-to-0 decay idx has been identified
                    else
                        if (data.final(y).smoothed_interval(idx_counter)~=0 && ~isnan(data.final(y).smoothed_interval(idx_counter)))
                            index = index + 1;
                        elseif (isnan(data.final(y).smoothed_interval(idx_counter)) && nan_counter == 0)
                            number_stim = number_stim + 1;
                            nan_counter = 1;
                        elseif (~isnan(data.final(y).smoothed_interval(idx_counter)) && nan_counter == 1)
                            nan_counter = 0;
                        end
                        idx_counter = idx_counter + 1;
                    end
                end

                if (idx_counter > size(data.final(y).smoothed_interval,1)*size(data.final(y).smoothed_interval,2) || data.final(y).max_peak_fraction(1,1) < 0.1)
                   rise_intersection_at_10_idx = NaN;
                else
                   rise_intersection_at_10_idx = data.offset_idx(1,y) + rise_to_peak * (data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;
                end

                % Intersection at 90% dF/F of peak amp when the signal is rising towards peak.
                rise_to_peak = []; counter1 = 0; idx_counter= 1; index = 0;
                nan_counter = 0; number_stim = 0;                           % keeps track of how many stimuli we have gone through
                tolerance = 0.015;                                          % Tolerance value: (%dF/F). If the value falls within this range
                adjacent_noise_tolerance = 2;                               % Prevents detection of the stimuli as peaks (which results in negative rise/decay time values)
                while (counter1 == 0 && idx_counter < size(data.final(y).smoothed_interval,1)*size(data.final(y).smoothed_interval,2))
                    if (data.final(y).smoothed_interval(idx_counter) > 0 && data.final(y).smoothed_interval(idx_counter)~=0 && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).max_peak_fraction(1,3))<=tolerance && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).smoothed_interval(idx_counter+data.smoothing_window/2))<adjacent_noise_tolerance)
                        rise_to_peak = index;
                        counter1 = 1; % Break the loop when signal-to-0 decay idx has been identified
                    elseif (data.final(y).smoothed_interval(idx_counter) > data.final(y).max_peak_fraction(1,3) && data.final(y).smoothed_interval(idx_counter)~=0 && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).max_peak_fraction(1,3))>tolerance && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).smoothed_interval(idx_counter+data.smoothing_window/2))<adjacent_noise_tolerance)
                        rise_to_peak = index;
                        counter1 = 1; % Break the loop when signal-to-0 decay idx has been identified
                    else
                        if (data.final(y).smoothed_interval(idx_counter)~=0 && ~isnan(data.final(y).smoothed_interval(idx_counter)))
                            index = index + 1;
                        elseif (isnan(data.final(y).smoothed_interval(idx_counter)) && nan_counter == 0)
                            number_stim = number_stim + 1;
                            nan_counter = 1;
                        elseif (~isnan(data.final(y).smoothed_interval(idx_counter)) && nan_counter == 1)
                            nan_counter = 0;
                        end
                        idx_counter = idx_counter + 1;
                    end
                end

                if (idx_counter > size(data.final(y).smoothed_interval,1)*size(data.final(y).smoothed_interval,2) || data.final(y).max_peak_fraction(1,3) < 0.9)
                   rise_intersection_at_90_idx = NaN;
                else
                rise_intersection_at_90_idx = data.offset_idx(1,y) + number_stim * data.pulse_width(y) * data.Fs +  rise_to_peak * (data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;
                end
                % Calculate 10%->90% Signal rise time
                data.final(y).ten_to_ninety_rise_time = (rise_intersection_at_90_idx - rise_intersection_at_10_idx)/data.Fs;

                % Decay Time Constant: (This is approximately when the signal reaches 37% from its max value)
                % Calculate the fall intersection at which we see signal baselining again after peak
                % (I.e. Intersection at 0% dF/F)
                decay_to_37 = []; counter2 = 0;
                tolerance = 0.015;                                                                          % Tolerance value: (%dF/F). If the value falls within this range 
                idx_counter = find(data.final(y).smoothed_peak(end)==data.final(y).smoothed_interval);      % Always start index counter from the last peak that was seen.
                adjacent_noise_tolerance = 2;                                                               % Prevents detection of the stimuli as peaks (which results in negative rise/decay time values)
                while (counter2 == 0 && idx_counter < size(data.final(y).smoothed_interval,1)*size(data.final(y).smoothed_interval,2))
                    if (data.final(y).smoothed_interval(idx_counter) > 0 && data.final(y).smoothed_interval(idx_counter)~=0 && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).max_peak_fraction(1,2))<=tolerance && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).smoothed_interval(idx_counter+data.smoothing_window/2))<adjacent_noise_tolerance)
                        decay_to_37 = index;
                        counter2 = 1; % Break the loop when signal-to-0 decay idx hsa been identified
                    elseif (data.final(y).smoothed_interval(idx_counter) < data.final(y).max_peak_fraction(1,2) && data.final(y).smoothed_interval(idx_counter)~=0 && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).max_peak_fraction(1,2))>tolerance && abs(data.final(y).smoothed_interval(idx_counter)-data.final(y).smoothed_interval(idx_counter+data.smoothing_window/2))<adjacent_noise_tolerance)    
                        decay_to_37 = index;
                        counter2 = 1; % Break the loop when signal-to-0 decay idx hsa been identified
                    else
                        if (data.final(y).smoothed_interval(idx_counter)~=0 && ~isnan(data.final(y).smoothed_interval(idx_counter)))
                            index = index + 1;
                        elseif (isnan(data.final(y).smoothed_interval(idx_counter)) && nan_counter == 0)
                            number_stim = number_stim + 1;
                            nan_counter = 1;
                        elseif (~isnan(data.final(y).smoothed_interval(idx_counter)) && nan_counter == 1)
                            nan_counter = 0;
                        end
                        idx_counter = idx_counter + 1;
                    end
                end

                if (idx_counter < size(data.final(y).smoothed_interval,1)*size(data.final(y).smoothed_interval,2))
                    decay_intersection_at_37_idx = data.offset_idx(find(~isnan(data.offset_idx(:,y)),1,'last'),y) + number_stim * data.pulse_width(y) * data.Fs + decay_to_37 * (data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;
                    % Calculate the Decay Time Constant (at 37% of peak signal)
                    if (~isempty(decay_to_37) & (decay_intersection_at_37_idx - data.final(y).smoothed_peak_locs(end))/data.Fs > 0 & data.final(y).max_peak_fraction(1,2) > 0.37)
                         data.final(y).decay_time_constant = (decay_intersection_at_37_idx - data.final(y).smoothed_peak_locs(end))/data.Fs;
                    else
                        data.final(y).decay_time_constant = NaN;
                    end
                else
                    data.final(y).decay_time_constant = NaN;
                end
                counter2 = 0; idx_counter = 0; index = 0;   % Reset counters

          end

        elseif (data.mouse == "fig1e")
            for y = 1: process_idx
                  for z = 1:data.num_stim(y)
                        if (data.num_stim(y) == z)
                        % PEAK
                        % AVERAGING/SMOOTHING
                        % Default: 1ms window size. 0.5ms overlay.
                        smoothed_interval = [];
                        % For the last stim, use the method of manually determining the interval of which to look at. (Default 2000ms interval)
                        extracted_interval = data.final(y).baselineFP_dFF(data.offset_idx(z,y):data.offset_idx(z,y)+data.search_range*(data.Fs/1000)-1);
                        for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(y).baselineFP_dFF(data.offset_idx(z,y):data.offset_idx(z,y)+data.search_range*(data.Fs/1000)-1))-data.smoothing_window
                            smoothed_interval(zz,1) = mean(extracted_interval(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                        end
                        % Align the Smoothed Interval to Real Intervals (i.e. the time from start of interval to time of 1st sampling) for precise indexing
                        data.final(y).smoothed_interval(:,z) = [smoothed_interval; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval),1)];
                        % Remove the stimulation itself from the data before determining the peak
                        % Later on, removed indicies (240 indicies) are manually added for determining peak location.
                        if (y == 2)
                            smoothed_interval(1:111,z) = 0;      % Masking stim width (11ms)
                            smoothed_interval(2451:2561,z) = 0;  % Masking stim width (11ms)
                            smoothed_interval(4951:5061,z) = 0;  % Masking stim width (11ms)        
                            smoothed_interval(7451:7561,z) = 0;  % Masking stim width (11ms)      
                        elseif (y == 3)
                            smoothed_interval(1:111,z) = 0;      % Masking stim width (11ms)
                            smoothed_interval(1201:1311,z) = 0;  % Masking stim width (11ms)
                            smoothed_interval(2451:2561,z) = 0;  % Masking stim width (11ms)
                            smoothed_interval(3701:3811,z) = 0;  % Masking stim width (11ms)
                        elseif (y == 4)
                            smoothed_interval(1:111,z) = 0;      % Masking stim width (11ms)
                            smoothed_interval(576:686,z) = 0;    % Masking stim width (11ms)
                            smoothed_interval(1201:1311,z) = 0;  % Masking stim width (11ms) 
                            smoothed_interval(1826:1936,z) = 0;  % Masking stim width (11ms)
                        end
                        smoothed_interval_work = [zeros(round(data.smoothing_window/2),1); smoothed_interval(smoothed_interval~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval(smoothed_interval~=0))-data.smoothing_window/2,1)];
                        data.final(y).smoothed_interval_idx(:,z) = [find(smoothed_interval~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval(smoothed_interval~=0)),1)];
                        % Remove all values that are the stimulus itself (Approx Stimuli Interval * 3)
                        smoothed_interval_work(1:ceil(data.max_point_within(1)*(data.Fs/1000)/(data.smoothing_window-data.smoothing_overlay))) = 0;
                        % Remove all zeros (i.e. non-sampled indicies from the vector)
                        smoothed_interval_work(smoothed_interval_work == 0) = NaN;
    
                        % For the last stim, use the same method of manually determining the interval of which to look at.
                        counter = 0; tolerance = 1;
                        while (counter == 0)
                            [data.final(y).smoothed_peak(:,z),data.final(y).smoothed_peak_idx(:,z)] = nanmax(smoothed_interval_work);
                            if (data.final(y).smoothed_peak_idx(:,z) > 1)
                                if abs(smoothed_interval_work(data.final(y).smoothed_peak_idx(:,z))-smoothed_interval_work(data.final(y).smoothed_peak_idx(:,z)-1)) > tolerance | find(find(data.final(y).smoothed_interval(:,z)==data.final(y).smoothed_peak(:,z))==data.final(y).smoothed_interval_idx(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2 > data.max_point_within(2)*(data.Fs/1000)
                                    smoothed_interval_work(data.final(y).smoothed_peak_idx(:,z)) = NaN;
                                else
                                    counter = 1;
                                end
                            else
                                counter = 1;
                            end
                        end 
                        counter = 0;
                        data.final(y).smoothed_peak_locs(:,z) = data.offset_idx(z,y) + find(find(data.final(y).smoothed_interval(:,z)==data.final(y).smoothed_peak(:,z))==data.final(y).smoothed_interval_idx(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;
    
                        % DIP
                        % AVERAGING/SMOOTHING
                        smoothed_interval_post = [];
                        extracted_interval_post = data.final(y).baselineFP_dFF(data.offset_idx(z,y):data.offset_idx(z,y)+data.search_range*(data.Fs/1000)-1);
                        for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(y).baselineFP_dFF(data.offset_idx(z,y):data.offset_idx(z,y)+data.search_range*(data.Fs/1000)-1))-data.smoothing_window
                            smoothed_interval_post(zz,1) = mean(extracted_interval_post(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                        end
                        data.final(y).smoothed_interval_post(:,z) = [smoothed_interval_post; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_post),1)];
                      
                        % Remove the stimulation itself from the data before determining the peak
                        % Later on, removed indicies (240 indicies) are manually added for determining peak location.
                        if (y == 2)
                            smoothed_interval_post(1:111,z) = 0;      % Masking stim width (11ms)
                            smoothed_interval_post(2451:2561,z) = 0;  % Masking stim width (11ms)
                            smoothed_interval_post(4951:5061,z) = 0;  % Masking stim width (11ms)        
                            smoothed_interval_post(7451:7561,z) = 0;  % Masking stim width (11ms)      
                        elseif (y == 3)
                            smoothed_interval_post(1:111,z) = 0;      % Masking stim width (11ms)
                            smoothed_interval_post(1201:1311,z) = 0;  % Masking stim width (11ms)
                            smoothed_interval_post(2451:2561,z) = 0;  % Masking stim width (11ms)
                            smoothed_interval_post(3701:3811,z) = 0;  % Masking stim width (11ms)
                        elseif (y == 4)
                            smoothed_interval_post(1:111,z) = 0;      % Masking stim width (11ms)
                            smoothed_interval_post(576:686,z) = 0;    % Masking stim width (11ms)
                            smoothed_interval_post(1201:1311,z) = 0;  % Masking stim width (11ms) 
                            smoothed_interval_post(1826:1936,z) = 0;  % Masking stim width (11ms)
                        end
                        
                        smoothed_interval_work_post = [zeros(round(data.smoothing_window/2),1); smoothed_interval_post(smoothed_interval_post~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_post(smoothed_interval_post~=0))-data.smoothing_window/2,1)];
                        data.final(y).smoothed_interval_idx_post(:,z) = [find(smoothed_interval_post~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_post(smoothed_interval_post~=0)),1)];
                        smoothed_interval_work_post(1:ceil(data.max_point_within(1)*(data.Fs/1000)/(data.smoothing_window-data.smoothing_overlay))) = 0;
                        smoothed_interval_work_post(smoothed_interval_work_post == 0) = NaN;
                        counter = 0; tolerance = 1;
                        while (counter == 0)
                            [data.final(y).smoothed_dip(:,z),data.final(y).smoothed_dip_idx(:,z)] = nanmin(smoothed_interval_work_post);
                            if abs(smoothed_interval_work_post(data.final(y).smoothed_dip_idx(:,z))-smoothed_interval_work_post(data.final(y).smoothed_dip_idx(:,z)-1)) > tolerance | find(find(data.final(y).smoothed_interval_post(:,z)==data.final(y).smoothed_dip(:,z))==data.final(y).smoothed_interval_idx_post(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2 > data.max_point_within(2)*20*(data.Fs/1000)
                                smoothed_interval_work_post(data.final(y).smoothed_dip_idx(:,z)) = NaN;
                            else
                                counter = 1;
                            end
                        end 
                        counter = 0;
                        data.final(y).smoothed_dip_locs(:,z) = data.offset_idx(z,y) + find(find(data.final(y).smoothed_interval_post(:,z)==data.final(y).smoothed_dip(:,z))==data.final(y).smoothed_interval_idx_post(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;
    
                        % DIP
                        % Smoothing an interval prior to the stimulation to acquire the baseline. (i.e. in comparision to the effect size of the pulse on the dip of the baseline fluorescence)
                        smoothed_interval_pre = [];
                        extracted_interval_pre = data.final(y).baselineFP_dFF(data.onset_idx(z,y)-20*(data.Fs/1000):data.onset_idx(z,y)-1);
                        for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(y).baselineFP_dFF(data.onset_idx(z,y)-20*(data.Fs/1000):data.onset_idx(z,y)-1))-data.smoothing_window
                            smoothed_interval_pre(zz,1) = mean(extracted_interval_pre(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                        end
                        data.final(y).smoothed_interval_pre(:,z) = [smoothed_interval_pre; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_pre),1)];
                        smoothed_interval_work_pre = [zeros(round(data.smoothing_window/2),1); smoothed_interval_pre(smoothed_interval_pre~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_pre(smoothed_interval_pre~=0))-data.smoothing_window/2,1)];
                        data.final(y).smoothed_interval_idx_pre(:,z) = [find(smoothed_interval_pre~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_pre(smoothed_interval_pre~=0)),1)];
                        smoothed_interval_work_pre(smoothed_interval_work_pre>20) = 0;
                        smoothed_interval_work_pre(smoothed_interval_work_pre<-20) = 0;
                        smoothed_interval_work_pre(smoothed_interval_work_pre == 0) = NaN;
                        data.final(y).smoothed_pre_stim(:,z) = nanmean(smoothed_interval_work_pre);                       
                        if (data.final(y).smoothed_pre_stim(:,z) - data.final(y).smoothed_dip(:,z) > 0)
                            data.final(y).dip_amplitude(:,z) = data.final(y).smoothed_pre_stim(:,z) - data.final(y).smoothed_dip(:,z);
                        else
                            data.final(y).dip_amplitude(:,z) = NaN; % If dip doesn't occur, than set its amplitude to 0.
                        end
                else
    
                        % PEAK
                        % AVERAGING/SMOOTHING
                        % Default: 1ms window size. 0.5ms overlay.
                        smoothed_interval = [];
                        % For the last stim, use the method of manually determining the interval of which to look at. (Default 2000ms interval)
                        extracted_interval = data.final(y).baselineFP_dFF(data.offset_idx(z,y):data.onset_idx(z+1,y)-1);
                        for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(y).baselineFP_dFF(data.offset_idx(z,y):data.onset_idx(z+1,y)-1))-data.smoothing_window
                            smoothed_interval(zz,1) = mean(extracted_interval(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                        end
                        % Align the Smoothed Interval to Real Intervals (i.e. the time from start of interval to time of 1st sampling) for precise indexing
                        data.final(y).smoothed_interval(:,z) = [smoothed_interval; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval),1)];
                        % Remove the stimulation itself from the data before determining the peak
                        % Later on, removed indicies (240 indicies) are manually added
                        % for determining peak location.
                        if (y == 2)
                            smoothed_interval(1:111,z) = 0;      % Masking stim width (11ms)
                            smoothed_interval(2451:2561,z) = 0;  % Masking stim width (11ms)
                            smoothed_interval(4951:5061,z) = 0;  % Masking stim width (11ms)        
                            smoothed_interval(7451:7561,z) = 0;  % Masking stim width (11ms)      
                        elseif (y == 3)
                            smoothed_interval(1:111,z) = 0;      % Masking stim width (11ms)
                            smoothed_interval(1201:1311,z) = 0;  % Masking stim width (11ms)
                            smoothed_interval(2451:2561,z) = 0;  % Masking stim width (11ms)
                            smoothed_interval(3701:3811,z) = 0;  % Masking stim width (11ms)
                        elseif (y == 4)
                            smoothed_interval(1:111,z) = 0;      % Masking stim width (11ms)
                            smoothed_interval(576:686,z) = 0;    % Masking stim width (11ms)
                            smoothed_interval(1201:1311,z) = 0;  % Masking stim width (11ms) 
                            smoothed_interval(1826:1936,z) = 0;  % Masking stim width (11ms)
                        end
    
                        smoothed_interval_work = [zeros(round(data.smoothing_window/2),1); smoothed_interval(smoothed_interval~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval(smoothed_interval~=0))-data.smoothing_window/2,1)];
                        data.final(y).smoothed_interval_idx(:,z) = [find(smoothed_interval~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval(smoothed_interval~=0)),1)];
                        % Remove all values that are the stimulus itself (Approx Stimuli Interval * 3)
                        smoothed_interval_work(1:ceil(data.max_point_within(1)*(data.Fs/1000)/(data.smoothing_window-data.smoothing_overlay))) = 0;
                        % Remove all zeros (i.e. non-sampled indicies from the vector)
                        smoothed_interval_work(smoothed_interval_work == 0) = NaN;
                        % For the last stim, use the same method of manually determining the interval of which to look at.
                        counter = 0; tolerance = 1;
                        while (counter == 0)
                            [data.final(y).smoothed_peak(:,z),data.final(y).smoothed_peak_idx(:,z)] = nanmax(smoothed_interval_work);
                            if (data.final(y).smoothed_peak_idx > 1)
                                if abs(smoothed_interval_work(data.final(y).smoothed_peak_idx(:,z))-smoothed_interval_work(data.final(y).smoothed_peak_idx(:,z)-1)) > tolerance | find(find(data.final(y).smoothed_interval(:,z)==data.final(y).smoothed_peak(:,z))==data.final(y).smoothed_interval_idx(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2 > data.max_point_within(2)*(data.Fs/1000)
                                    smoothed_interval_work(data.final(y).smoothed_peak_idx(:,z)) = NaN;
                                else
                                    counter = 1;
                                end
                            else
                                counter = 1;
                            end
                        end 
                        counter = 0;
                        data.final(y).smoothed_peak_locs(:,z) = data.offset_idx(z,y) + find(find(data.final(y).smoothed_interval(:,z)==data.final(y).smoothed_peak(:,z))==data.final(y).smoothed_interval_idx(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;
    
                       % DIP
                       % AVERAGING/SMOOTHING
                       smoothed_interval_post = [];
                       extracted_interval_post = data.final(y).baselineFP_dFF(data.offset_idx(z,y):data.onset_idx(z+1,y));
                       for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(1).baselineFP_dFF(data.offset_idx(z,y):data.onset_idx(z+1,y)))-data.smoothing_window
                           smoothed_interval_post(zz,1) = mean(extracted_interval_post(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                       end
                        data.final(y).smoothed_interval_post(:,z) = [smoothed_interval_post; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_post),1)];
                        smoothed_interval_work_post = [zeros(round(data.smoothing_window/2),1); smoothed_interval_post(smoothed_interval_post~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_post(smoothed_interval_post~=0))-data.smoothing_window/2,1)];
                        data.final(y).smoothed_interval_idx_post(:,z) = [find(smoothed_interval_post~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_post(smoothed_interval_post~=0)),1)];
                        smoothed_interval_work_post(1:ceil(data.max_point_within(1)*(data.Fs/1000)/(data.smoothing_window-data.smoothing_overlay))) = 0;
                        smoothed_interval_work_post(smoothed_interval_work_post == 0) = NaN;
                        counter = 0; tolerance = 1;
                        while (counter == 0)
                            [data.final(y).smoothed_dip(:,z),data.final(y).smoothed_dip_idx(:,z)] = nanmin(smoothed_interval_work_post);
                            if abs(smoothed_interval_work_post(data.final(y).smoothed_dip_idx(:,z))-smoothed_interval_work_post(data.final(y).smoothed_dip_idx(:,z)-1)) > tolerance | find(find(data.final(y).smoothed_interval_post(:,z)==data.final(y).smoothed_dip(:,z))==data.final(y).smoothed_interval_idx_post(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2 > data.max_point_within(2)*(data.Fs/1000)
                                smoothed_interval_work_post(data.final(y).smoothed_dip_idx(:,z)) = NaN;
                            else
                                counter = 1;
                            end
                        end 
                        counter = 0;
                        data.final(y).smoothed_dip_locs(:,z) = data.offset_idx(z,y) + find(find(data.final(y).smoothed_interval_post(:,z)==data.final(y).smoothed_dip(:,z))==data.final(y).smoothed_interval_idx_post(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;
                        
                        % Smoothing an interval prior to the stimulation to acquire the baseline. (i.e. in comparision to the effect size of the pulse on the dip of the baseline fluorescence)
                        smoothed_interval_pre = [];
                        extracted_interval_pre = data.final(y).baselineFP_dFF(data.onset_idx(z,y)-20*(data.Fs/1000):data.onset_idx(z,y)-1);
                        for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(y).baselineFP_dFF(data.onset_idx(z,y)-20*(data.Fs/1000):data.onset_idx(z,y)-1))-data.smoothing_window
                            smoothed_interval_pre(zz,1) = mean(extracted_interval_pre(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                        end
                        data.final(y).smoothed_interval_pre(:,z) = [smoothed_interval_pre; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_pre),1)];
                        smoothed_interval_work_pre = [zeros(round(data.smoothing_window/2),1); smoothed_interval_pre(smoothed_interval_pre~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_pre(smoothed_interval_pre~=0))-data.smoothing_window/2,1)];
                        data.final(y).smoothed_interval_idx_pre(:,z) = [find(smoothed_interval_pre~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_pre(smoothed_interval_pre~=0)),1)];
                        smoothed_interval_work_pre(smoothed_interval_work_pre>20) = 0;
                        smoothed_interval_work_pre(smoothed_interval_work_pre<-20) = 0;
                        smoothed_interval_work_pre(smoothed_interval_work_pre == 0) = NaN;
                        data.final(y).smoothed_pre_stim(:,z) = nanmean(smoothed_interval_work_pre);
                        if (data.final(y).smoothed_pre_stim(:,z) - data.final(y).smoothed_dip(:,z) > 0)
                            data.final(y).dip_amplitude(:,z) = data.final(y).smoothed_pre_stim(:,z) - data.final(y).smoothed_dip(:,z);
                        else
                            data.final(y).dip_amplitude(:,z) = NaN; % If dip doesn't occur, than set its amplitude to NaN.
                        end
                    end
                  
                  % Calculate the max dip amplitude (relative to each stim) across all
                  data.final(y).max_dip_amplitude = nanmax(data.final(y).dip_amplitude);
                  if data.final(y).max_dip_amplitude < 0
                     data.final(y).max_dip_amplitude = NaN; % NaN if no dips occur
                  end
                  % Calculate the mean dip amplitude per each condition
                  data.final(y).mean_dip_amplitude = nanmean(data.final(y).dip_amplitude);
                  % Calculate the total (i.e. accumulated dip) across all 
                  data.final(y).total_dip = data.final(y).smoothed_pre_stim(1)-min(data.final(y).smoothed_dip);
                  if data.final(y).total_dip < 0
                     data.final(y).total_dip = NaN; % NaN if no dips occur
                  end
              end
        end
  

        elseif (data.mouse == "fig1f")
            for y = 1:process_idx
                for z = 1:data.num_stim(y)
                    if (data.num_stim(y) == z)
                        % PEAK
                        % AVERAGING/SMOOTHING
                        % Default: 1ms window size. 0.5ms overlay.
                        smoothed_interval = [];
                        % For the last stim, use the method of manually determining the interval of which to look at. (Default 2000ms interval)
                        extracted_interval = data.final(y).baselineFP_dFF(data.offset_idx(z,y):data.offset_idx(z,y)+data.search_range*(data.Fs/1000)-1);
                        for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(y).baselineFP_dFF(data.offset_idx(z,y):data.offset_idx(z,y)+data.search_range*(data.Fs/1000)-1))-data.smoothing_window
                            smoothed_interval(zz,1) = mean(extracted_interval(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                        end
                        % Align the Smoothed Interval to Real Intervals (i.e. the time from start of interval to time of 1st sampling) for precise indexing
                        data.final(y).smoothed_interval(:,z) = [smoothed_interval; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval),1)];
                        smoothed_interval_work = [zeros(round(data.smoothing_window/2),1); smoothed_interval(smoothed_interval~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval(smoothed_interval~=0))-data.smoothing_window/2,1)];
                        data.final(y).smoothed_interval_idx(:,z) = [find(smoothed_interval~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval(smoothed_interval~=0)),1)];
                        % Remove all values that are the stimulus itself (Approx Stimuli Interval * 3)
                        smoothed_interval_work(1:ceil(data.max_point_within(1)*(data.Fs/1000)/(data.smoothing_window-data.smoothing_overlay))) = 0;
                        % Remove all zeros (i.e. non-sampled indicies from the vector)
                        smoothed_interval_work(smoothed_interval_work == 0) = NaN;

                        % For the last stim, use the same method of manually determining the interval of which to look at.
                        % Counter prevents detecting smoothed artifacts as peaks
                        counter = 0; tolerance = 1;
                        while (counter == 0)
                                [data.final(y).smoothed_peak(:,z),data.final(y).smoothed_peak_idx(:,z)] = nanmax(smoothed_interval_work);
                            if (data.final(y).smoothed_peak_idx(:,z) > 1)
                                if abs(smoothed_interval_work(data.final(y).smoothed_peak_idx(:,z))-smoothed_interval_work(data.final(y).smoothed_peak_idx(:,z)-1)) > tolerance | find(find(data.final(y).smoothed_interval(:,z)==data.final(y).smoothed_peak(:,z))==data.final(y).smoothed_interval_idx(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2 > data.max_point_within(2)*(data.Fs/1000)
                                    smoothed_interval_work(data.final(y).smoothed_peak_idx(:,z)) = NaN;
                                else
                                    counter = 1;
                                end
                            else
                                counter = 1;
                            end
                        end 
                        counter = 0;
                        data.final(y).smoothed_peak_locs(:,z) = data.offset_idx(z,y) + find(find(data.final(y).smoothed_interval(:,z)==data.final(y).smoothed_peak(:,z))==data.final(y).smoothed_interval_idx(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;
                    
                        % DIP
                        % AVERAGING/SMOOTHING
                        smoothed_interval_post = [];
                        extracted_interval_post = data.final(y).baselineFP_dFF(data.offset_idx(z,y):data.offset_idx(z,y)+data.search_range*(data.Fs/1000)-1);
                        for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(y).baselineFP_dFF(data.offset_idx(z,y):data.offset_idx(z,y)+data.search_range*(data.Fs/1000)-1))-data.smoothing_window
                            smoothed_interval_post(zz,1) = mean(extracted_interval_post(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                        end
                        data.final(y).smoothed_interval_post(:,z) = [smoothed_interval_post; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_post),1)];
                        smoothed_interval_work_post = [zeros(round(data.smoothing_window/2),1); smoothed_interval_post(smoothed_interval_post~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_post(smoothed_interval_post~=0))-data.smoothing_window/2,1)];
                        data.final(y).smoothed_interval_idx_post(:,z) = [find(smoothed_interval_post~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_post(smoothed_interval_post~=0)),1)];
                        smoothed_interval_work_post(1:ceil(data.max_point_within(1)*(data.Fs/1000)/(data.smoothing_window-data.smoothing_overlay))) = 0;
                        smoothed_interval_work_post(smoothed_interval_work_post == 0) = NaN;
                        counter = 0; tolerance = 1;
                        while (counter == 0)
                            [data.final(y).smoothed_dip(:,z),data.final(y).smoothed_dip_idx(:,z)] = nanmin(smoothed_interval_work_post);
                            if abs(smoothed_interval_work_post(data.final(y).smoothed_dip_idx(:,z))-smoothed_interval_work_post(data.final(y).smoothed_dip_idx(:,z)-1)) > tolerance | find(find(data.final(y).smoothed_interval_post(:,z)==data.final(y).smoothed_dip(:,z))==data.final(y).smoothed_interval_idx_post(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2 > data.max_point_within(2)*(data.Fs/1000)
                                smoothed_interval_work_post(data.final(y).smoothed_dip_idx(:,z)) = NaN;
                            else
                                counter = 1;
                            end
                        end 
                        counter = 0;
                        data.final(y).smoothed_dip_locs(:,z) = data.offset_idx(z,y) + find(find(data.final(y).smoothed_interval_post(:,z)==data.final(y).smoothed_dip(:,z))==data.final(y).smoothed_interval_idx_post(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;
                              
                        % DIP
                        % Smoothing an interval prior to the stimulation to acquire the baseline. (i.e. in comparision
                        % to the effect size of the pulse on the dip of the baseline fluorescence)
                        smoothed_interval_pre = [];
                        extracted_interval_pre = data.final(y).baselineFP_dFF(data.onset_idx(z,y)-20*(data.Fs/1000):data.onset_idx(z,y)-1);
                        for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(y).baselineFP_dFF(data.onset_idx(z,y)-20*(data.Fs/1000):data.onset_idx(z,y)-1))-data.smoothing_window
                            smoothed_interval_pre(zz,1) = mean(extracted_interval_pre(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                        end
                        data.final(y).smoothed_interval_pre(:,z) = [smoothed_interval_pre; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_pre),1)];
                        smoothed_interval_work_pre = [zeros(round(data.smoothing_window/2),1); smoothed_interval_pre(smoothed_interval_pre~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_pre(smoothed_interval_pre~=0))-data.smoothing_window/2,1)];
                        data.final(y).smoothed_interval_idx_pre(:,z) = [find(smoothed_interval_pre~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_pre(smoothed_interval_pre~=0)),1)];
                        smoothed_interval_work_pre(smoothed_interval_work_pre>20) = 0;
                        smoothed_interval_work_pre(smoothed_interval_work_pre<-20) = 0;
                        smoothed_interval_work_pre(smoothed_interval_work_pre == 0) = NaN;
                        data.final(y).smoothed_pre_stim(:,z) = nanmean(smoothed_interval_work_pre);
                        if (data.final(y).smoothed_pre_stim(:,z) - data.final(y).smoothed_dip(:,z) > 0)
                            data.final(y).dip_amplitude(:,z) = data.final(y).smoothed_pre_stim(:,z) - data.final(y).smoothed_dip(:,z);
                        else
                            data.final(y).dip_amplitude(:,z) = NaN; % If dip doesn't occur, than set its amplitude to 0.
                        end
                    else 

                       % PEAK
                       % AVERAGING/SMOOTHING
                       % Default: 1ms window size. 0.5ms overlay.
                       smoothed_interval = [];
                       extracted_interval = data.final(y).baselineFP_dFF(data.offset_idx(z,y):data.offset_idx(z,y)+data.search_range*(data.Fs/1000)-1);
                       for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(y).baselineFP_dFF(data.offset_idx(z,y):data.offset_idx(z,y)+data.search_range*(data.Fs/1000)-1))-data.smoothing_window
                           smoothed_interval(zz,1) = mean(extracted_interval(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                       end

                       % Align the Smoothed Interval to Real Intervals (i.e. the time from start of interval to time of 1st sampling) for precise indexing
                       data.final(y).smoothed_interval(:,z) = [smoothed_interval; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval),1)];
                       smoothed_interval_work = [zeros(round(data.smoothing_window/2),1); smoothed_interval(smoothed_interval~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval(smoothed_interval~=0))-data.smoothing_window/2,1)];
                       data.final(y).smoothed_interval_idx(:,z) = [find(smoothed_interval~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval(smoothed_interval~=0)),1)];
                       % Remove all values that are the stimulus itself (Approx Stimuli Interval * 3)
                       smoothed_interval_work(1:ceil(data.max_point_within(1)*(data.Fs/1000)/(data.smoothing_window-data.smoothing_overlay))) = 0;
                       % Remove all zeros (i.e. non-sampled indicies from the vector)
                       smoothed_interval_work(smoothed_interval_work == 0) = NaN;

                       % For the last stim, use the same method of manually determining the interval of which to look at.
                       % Counter prevents detecting smoothed artifacts as peaks
                       counter = 0; tolerance = 1;
                       while (counter == 0)
                           [data.final(y).smoothed_peak(:,z),data.final(y).smoothed_peak_idx(:,z)] = nanmax(smoothed_interval_work);
                           if (data.final(y).smoothed_peak_idx > 1)
                               if abs(smoothed_interval_work(data.final(y).smoothed_peak_idx(:,z))-smoothed_interval_work(data.final(y).smoothed_peak_idx(:,z)-1)) > tolerance | find(find(data.final(y).smoothed_interval(:,z)==data.final(y).smoothed_peak(:,z))==data.final(y).smoothed_interval_idx(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2 > data.max_point_within(2)*(data.Fs/1000)
                                   smoothed_interval_work(data.final(y).smoothed_peak_idx(:,z)) = NaN;
                               else
                                   counter = 1;
                               end
                           else
                               counter = 1;
                           end
                       end
                       counter = 0;
                       data.final(y).smoothed_peak_locs(:,z) = data.offset_idx(z,y) + find(find(data.final(y).smoothed_interval(:,z)==data.final(y).smoothed_peak(:,z))==data.final(y).smoothed_interval_idx(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;

                       % DIP
                       % AVERAGING/SMOOTHING
                       smoothed_interval_post = [];
                       extracted_interval_post = data.final(y).baselineFP_dFF(data.offset_idx(z,y):data.offset_idx(z,y)+data.search_range*(data.Fs/1000));
                       for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(1).baselineFP_dFF(data.offset_idx(z,y):data.offset_idx(z,y)+data.search_range*(data.Fs/1000)))-data.smoothing_window
                           smoothed_interval_post(zz,1) = mean(extracted_interval_post(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                       end
                       data.final(y).smoothed_interval_post(:,z) = [smoothed_interval_post; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_post),1)];
                       smoothed_interval_work_post = [zeros(round(data.smoothing_window/2),1); smoothed_interval_post(smoothed_interval_post~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_post(smoothed_interval_post~=0))-data.smoothing_window/2,1)];
                       data.final(y).smoothed_interval_idx_post(:,z) = [find(smoothed_interval_post~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_post(smoothed_interval_post~=0)),1)];
                       smoothed_interval_work_post(1:ceil(data.max_point_within(1)*(data.Fs/1000)/(data.smoothing_window-data.smoothing_overlay))) = 0;
                       smoothed_interval_work_post(smoothed_interval_work_post == 0) = NaN;
                       counter = 0; tolerance = 1;
                       while (counter == 0)
                           [data.final(y).smoothed_dip(:,z),data.final(y).smoothed_dip_idx(:,z)] = nanmin(smoothed_interval_work_post);
                           if abs(smoothed_interval_work_post(data.final(y).smoothed_dip_idx(:,z))-smoothed_interval_work_post(data.final(y).smoothed_dip_idx(:,z)-1)) > tolerance | find(find(data.final(y).smoothed_interval_post(:,z)==data.final(y).smoothed_dip(:,z))==data.final(y).smoothed_interval_idx_post(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2 > data.max_point_within(2)*(data.Fs/1000)
                               smoothed_interval_work_post(data.final(y).smoothed_dip_idx(:,z)) = NaN;
                           else
                               counter = 1;
                           end
                       end 
                       counter = 0;
                       data.final(y).smoothed_dip_locs(:,z) = data.offset_idx(z,y) + find(find(data.final(y).smoothed_interval_post(:,z)==data.final(y).smoothed_dip(:,z))==data.final(y).smoothed_interval_idx_post(:,z))*(data.smoothing_window-data.smoothing_overlay)+data.smoothing_window/2;
                       
                       % DIP
                       % Smoothing an interval prior to the stimulation to acquire the baseline. (i.e. in comparision
                       % to the effect size of the pulse on the dip of the baseline fluorescence)
                       smoothed_interval_pre = [];
                       extracted_interval_pre = data.final(y).baselineFP_dFF(data.onset_idx(z,y)-20*(data.Fs/1000):data.onset_idx(z,y)-1);
                       for zz = 1:data.smoothing_window-data.smoothing_overlay:length(data.final(y).baselineFP_dFF(data.onset_idx(z,y)-20*(data.Fs/1000):data.onset_idx(z,y)-1))-data.smoothing_window
                           smoothed_interval_pre(zz,1) = mean(extracted_interval_pre(zz+round(data.smoothing_window/2)-round(data.smoothing_window/2):zz+round(data.smoothing_window/2)+round(data.smoothing_window/2)));
                       end
                       data.final(y).smoothed_interval_pre(:,z) = [smoothed_interval_pre; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_pre),1)];
                       smoothed_interval_work_pre = [zeros(round(data.smoothing_window/2),1); smoothed_interval_pre(smoothed_interval_pre~=0); nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_pre(smoothed_interval_pre~=0))-data.smoothing_window/2,1)];
                       data.final(y).smoothed_interval_idx_pre(:,z) = [find(smoothed_interval_pre~=0)+data.smoothing_window/2; nan(data.search_range*(data.Fs/1000)-length(smoothed_interval_pre(smoothed_interval_pre~=0)),1)];
                       smoothed_interval_work_pre(smoothed_interval_work_pre>20) = 0;
                       smoothed_interval_work_pre(smoothed_interval_work_pre<-20) = 0;
                       smoothed_interval_work_pre(smoothed_interval_work_pre == 0) = NaN;
                       data.final(y).smoothed_pre_stim(:,z) = nanmean(smoothed_interval_work_pre);
                       if (data.final(y).smoothed_pre_stim(:,z) - data.final(y).smoothed_dip(:,z) > 0)
                           data.final(y).dip_amplitude(:,z) = data.final(y).smoothed_pre_stim(:,z) - data.final(y).smoothed_dip(:,z);
                       else
                           data.final(y).dip_amplitude(:,z) = NaN; % If dip doesn't occur, than set its amplitude to 0.
                       end
                    end
                end
            end
        end
    % Save the data matrix
    save(fullfile(FPpath,FPfiles{x}),'data','-v7.3');
end
end
