function [autocorr_ACh, autocorr_DA, corr_AChACh, corr_DADA, corr_cell_samehemi, lags, shuffAChACh, shuffDADA, shuff_achda] = JT_autocorrFP(beh, varargin)
%Cross-correlation between ACh and DA photometry signals acquired during
%dual color photometry recordings
%
%   [corr_achda, lags] = AK_corrFP(beh)
%   [corr_achda, lags] = AK_corrFP(beh, win)
%   [corr_achda, lags, shuff_achda, min_val, min_lag] = AK_corrFP(beh)
%
%   Description: This function is for running cross-correlation analysis on
%   two continuous photometry signals (ACh, DA) after isolating time
%   periods when animal is immobile (no movement, no reward) using the
%   MATLAB function 'xcorr'
%
%   INPUTS
%   'beh' - structure with photometry and behavioral data for multiple
%   recordings, should include beh(x).rec as [an,'-',day]
%   'win'(optional) - window to restrict analysis to, in seconds
%   
%   OUTPUTS
%   'corr_achda' - cross-correlation between ACh and DA signals, normalized
%   using coeff scaling to approximate Pearson's coefficient, r.
%   'lags' - vector of lag indices from +/- 5 seconds, in samples
%   'shuff_achda' - 5th, 50th, 95th percentiles for cross-correlation run
%   on shuffled photometry signals, across behavioral states
%       shuff_achda{1,1} is 5th percentile for immobility, and so on
%       shuff_achda{3,2} is 50th percentile for reward, and so on
%       shuff_achda{2,3} is 95th percentile for locomotion, and so on
%
%   Author: Anya Krok, December 2021

%% INPUTS
if nargin == 2; win = varargin{1}; end % Window for analysis, in seconds
%y = [2 1]; % Photometry signal to use as reference is the one listed first
% e.g. if y = [1 2] then the signal in beh(x).FP{1} will be used as
% reference signal, while if y = [2 1] then the signal in beh(x).FP{2} will
% be used as reference signal instead
% Default from Krok 2022 is to use y = [2 1] so that rDA1m photometry
% signal is the reference signal
winCorr = 2; % Window, in seconds, to run cross-correlation on

analysis_type = menu('2 or 4 channels?','2','4');
channel_type = menu('DA/ACh or DA/DA or ACh/ACh ?','DA/ACh','same/same');
behavior_type = menu('Behavior Type?','Wheel', 'Freely Behaving?');
%animal_type = menu('Tony's or mine?','Tony','Mine');


%% OUTPUTS
% mat = struct; % temporary output structure

if (analysis_type==1)
    corr_cell_acrosshemi = cell(3,4); % temporary output cell array
    for a = 1:3; for b = 1:4; corr_cell_acrosshemi{a,b} = nan((winCorr*2*50)+1,length(beh)); end; end % fill cell array with nan's
else
    corr_cell_samehemi = cell(3,4); % temporary output cell array
    for a = 1:3; for b = 1:4; corr_cell_samehemi{a,b} = nan((winCorr*2*50)+1,length(beh)); end; end % fill cell array with nan's
end


%% RUN ANALYSIS ON ALL RECORDINGS
h = waitbar(0, 'cross-correlation and auto-correlation');
idxStates = extractBehavioralStates(beh);
nStates = 3;
counter=1;
for x = 1:length(beh) % iterate over all recordings
    
    %% extract signals

    % (NOT IN USE) For bilateral (3 or 4 channels), corr between DA and ACh:
    %for xx = 1:length(beh(x).FP)/2
    % fp_mat = [];
    % fp_mat(:,2*xx-1) = beh(x).FP{y(2*xx-1)}; Fs = beh(x).Fs; % extract photometry signal from structure, which will be used as reference
    % fp_mat(:,2*xx-1) = fp_mat(:,2*xx-1) - nanmean(fp_mat(:,2*xx-1)); % subtract baseline (mean of entire photometry signal) from fp
    % fp_mat(:,2*xx) = beh(x).FP{y(2*xx)}; % extract photometry signal from structure
    % fp_mat(:,2*xx) = fp_mat(:,2*xx) - nanmean(fp_mat(:,2*xx)); % subtract baseline (mean of entire photometry signal) from fp
  
    % Normal DA vs ACh 2-channel xcorr
    % DA vs DA (One Channel only)
     %fp_mat = [];
     %y = [2 1];
     %fp_mat(:,1) = beh(x).FP{y(1)}; Fs = beh(x).Fs; % extract photometry signal from structure, which will be used as reference
     %fp_mat(:,1) = fp_mat(:,1) - nanmean(fp_mat(:,1)); % subtract baseline (mean of entire photometry signal) from fp
     %fp_mat(:,2) = beh(x).FP{y(2)}; % extract photometry signal from structure
     %fp_mat(:,2) = fp_mat(:,2) - nanmean(fp_mat(:,2)); % subtract baseline (mean of entire photometry signal) from fp


    %% Bilateral DA vs ACh 2-channel xcorr

    % Corr across hemisphere (i.e. DA vs DA and ACh vs ACh)
   for xx = 1:length(beh(x).FP)/2
    if (analysis_type==2)
     yy = [xx+2 xx];
     fp_mat_acrosshemi = [];
     fp_mat_acrosshemi(:,1) = beh(x).FP{yy(1)}; Fs = beh(x).Fs; % extract photometry signal from structure, which will be used as reference
     fp_mat_acrosshemi(:,1) = fp_mat_acrosshemi(:,1) - nanmean(fp_mat_acrosshemi(:,1)); % subtract baseline (mean of entire photometry signal) from fp
     fp_mat_acrosshemi(:,2) = beh(x).FP{yy(2)}; % extract photometry signal from structure
     fp_mat_acrosshemi(:,2) = fp_mat_acrosshemi(:,2) - nanmean(fp_mat_acrosshemi(:,2)); % subtract baseline (mean of entire photometry signal) from fp
    end
    % DA vs ACh (One Channel only)
     fp_mat_samehemi = [];
     y = [2*xx 2*xx-1];
     fp_mat_samehemi(:,1) = beh(x).FP{y(1)}; Fs = beh(x).Fs; % extract photometry signal from structure, which will be used as reference
     fp_mat_samehemi(:,1) = fp_mat_samehemi(:,1) - nanmean(fp_mat_samehemi(:,1)); % subtract baseline (mean of entire photometry signal) from fp
     fp_mat_samehemi(:,2) = beh(x).FP{y(2)}; % extract photometry signal from structure
     fp_mat_samehemi(:,2) = fp_mat_samehemi(:,2) - nanmean(fp_mat_samehemi(:,2)); % subtract baseline (mean of entire photometry signal) from fp

    
    % 230608 UPDATE: If taking the 1st derivative of the raw signals
    %fp_mat_diff=[];
    %fp_mat_diff(:,1)=diff(fp_mat(:,1));
    %fp_mat_diff(:,2)=diff(fp_mat(:,2));
    %fp_mat=fp_mat_diff;
    
    
    %% adjust indices to retain if within specified window
%     if nargin == 2
%         for z = 1:nStates
%             idxTmp = idxStates{x,z};
%             idxTmp = idxTmp(idxTmp > win(1)*Fs & idxTmp < win(2)*Fs); % retain only indices that are within specified window
%             if ~isempty(idxTmp)
%                 needL = 200.*Fs;
%                 if length(idxTmp) < needL % ensure that have at least 200 seconds of signal
%                     idxTmp = repmat(idxTmp, [ceil(needL/length(idxTmp)) 1]); % duplicate indices to lengthen signal for processing
%                 end
%             end
%             idxStates{x,z} = idxTmp; % re-insert into output structure
%         end
%     end
    
    %%
    fp_sub_acrosshemi = [];
    fp_sub_samehemi=[];

        
        if (behavior_type==1)
          for z = 1:nStates
                if length(idxStates{x,z})< 2; continue; end % Only calculate states that have more than 2 indicies
                
                if (analysis_type==2)
                    % Cross-Correlation
                    fp_sub_acrosshemi = fp_mat_acrosshemi(idxStates{x,z},:); % signal
                    [corr_tmp_acrosshemi, lags] = xcorr(fp_sub_acrosshemi(:,1), fp_sub_acrosshemi(:,2), winCorr*Fs, 'coeff'); % cross-correlation
                   
                    fp_sub_new_acrosshemi = fp_sub_acrosshemi(:,2);
                    tmp_shuff_acrosshemi = []; 
                    for s = 1:50
                        fp_sub_new_acrosshemi = circshift(fp_sub_new_acrosshemi, Fs);
                        tmp_shuff_acrosshemi(:,s) = xcorr(fp_sub_acrosshemi(:,1), fp_sub_new_acrosshemi, winCorr*Fs, 'coeff');
                    end
        
                    % Corr across hemisphere (i.e. DA vs DA and ACh vs ACh)
                    corr_cell_acrosshemi{z,1}(:,counter) = corr_tmp_acrosshemi;                  % cross-correlation
                    corr_cell_acrosshemi{z,2}(:,counter) = prctile(tmp_shuff_acrosshemi, 5, 2);  % shuffle 5th percentile
                    corr_cell_acrosshemi{z,3}(:,counter) = prctile(tmp_shuff_acrosshemi, 50, 2); % shuffle 50th percentile
                    corr_cell_acrosshemi{z,4}(:,counter) = prctile(tmp_shuff_acrosshemi, 95, 2); % shuffle 95th percentile
                end
        
                % Cross-Correlation
                fp_sub_samehemi = fp_mat_samehemi(idxStates{x,z},:); % signal
                [corr_tmp_samehemi, lags] = xcorr(fp_sub_samehemi(:,1), fp_sub_samehemi(:,2), winCorr*Fs, 'coeff'); % cross-correlation
                
                fp_sub_new_samehemi = fp_sub_samehemi(:,2);
                tmp_shuff_samehemi = []; 
                for s = 1:50
                    fp_sub_new_samehemi = circshift(fp_sub_new_samehemi, Fs);
                    tmp_shuff_samehemi(:,s) = xcorr(fp_sub_samehemi(:,1), fp_sub_new_samehemi, winCorr*Fs, 'coeff');
                end
        
                % Normal DA vs ACh 2-channel xcorr
                corr_cell_samehemi{z,1}(:,counter) = corr_tmp_samehemi;                     % cross-correlation
                corr_cell_samehemi{z,2}(:,counter) = prctile(tmp_shuff_samehemi, 5, 2);     % shuffle 5th percentile
                corr_cell_samehemi{z,3}(:,counter) = prctile(tmp_shuff_samehemi, 50, 2);    % shuffle 50th percentile
                corr_cell_samehemi{z,4}(:,counter) = prctile(tmp_shuff_samehemi, 95, 2);    % shuffle 95th percentile
        
                if (analysis_type==2)
                    % Auto-correlation
                    [autocorr_tmp1,lags_tmp1]=autocorr(fp_sub_acrosshemi(:,1),'NumLags',winCorr*Fs);
                    [autocorr_tmp2,lags_tmp2]=autocorr(fp_sub_acrosshemi(:,2),'NumLags',winCorr*Fs);
                else
                    % Auto-correlation
                    [autocorr_tmp1,lags_tmp1]=autocorr(fp_sub_samehemi(:,1),'NumLags',winCorr*Fs);
                    [autocorr_tmp2,lags_tmp2]=autocorr(fp_sub_samehemi(:,2),'NumLags',winCorr*Fs);
                end
        
                % Autocorr for DA-ACh pairs channels (saves ACh from autocorr_tmp1 and DA from autocorr_tmp2)
                autocorr_cell{z,1}(:,1)=lags_tmp1;
                autocorr_cell{z,1}(:,2*counter)=autocorr_tmp1;
                autocorr_cell{z,1}(:,2*counter+1)=autocorr_tmp2;
          end
        else
            if (analysis_type==2)
                % Cross-Correlation
                fp_sub_acrosshemi = fp_mat_acrosshemi; % signal
                [corr_tmp_acrosshemi, lags] = xcorr(fp_sub_acrosshemi(:,1), fp_sub_acrosshemi(:,2), winCorr*Fs, 'coeff'); % cross-correlation
               
                fp_sub_new_acrosshemi = fp_sub_acrosshemi(:,2);
                tmp_shuff_acrosshemi = []; 
                for s = 1:50
                    fp_sub_new_acrosshemi = circshift(fp_sub_new_acrosshemi, Fs);
                    tmp_shuff_acrosshemi(:,s) = xcorr(fp_sub_acrosshemi(:,1), fp_sub_new_acrosshemi, winCorr*Fs, 'coeff');
                end
    
                % Corr across hemisphere (i.e. DA vs DA and ACh vs ACh)
                corr_cell_acrosshemi{1,1}(:,counter) = corr_tmp_acrosshemi;                  % cross-correlation
                corr_cell_acrosshemi{1,2}(:,counter) = prctile(tmp_shuff_acrosshemi, 5, 2);  % shuffle 5th percentile
                corr_cell_acrosshemi{1,3}(:,counter) = prctile(tmp_shuff_acrosshemi, 50, 2); % shuffle 50th percentile
                corr_cell_acrosshemi{1,4}(:,counter) = prctile(tmp_shuff_acrosshemi, 95, 2); % shuffle 95th percentile
            end
    
            % Cross-Correlation
            fp_sub_samehemi = fp_mat_samehemi(idxStates{x,1},:); % signal
            [corr_tmp_samehemi, lags] = xcorr(fp_sub_samehemi(:,1), fp_sub_samehemi(:,2), winCorr*Fs, 'coeff'); % cross-correlation
            
            fp_sub_new_samehemi = fp_sub_samehemi(:,2);
            tmp_shuff_samehemi = []; 
            for s = 1:50
                fp_sub_new_samehemi = circshift(fp_sub_new_samehemi, Fs);
                tmp_shuff_samehemi(:,s) = xcorr(fp_sub_samehemi(:,1), fp_sub_new_samehemi, winCorr*Fs, 'coeff');
            end
    
            % Normal DA vs ACh 2-channel xcorr
            corr_cell_samehemi{1,1}(:,counter) = corr_tmp_samehemi;                     % cross-correlation
            corr_cell_samehemi{1,2}(:,counter) = prctile(tmp_shuff_samehemi, 5, 2);     % shuffle 5th percentile
            corr_cell_samehemi{1,3}(:,counter) = prctile(tmp_shuff_samehemi, 50, 2);    % shuffle 50th percentile
            corr_cell_samehemi{1,4}(:,counter) = prctile(tmp_shuff_samehemi, 95, 2);    % shuffle 95th percentile
    
            if (analysis_type==2)
                % Auto-correlation
                [autocorr_tmp1,lags_tmp1]=autocorr(fp_sub_acrosshemi(:,1),'NumLags',winCorr*Fs);
                [autocorr_tmp2,lags_tmp2]=autocorr(fp_sub_acrosshemi(:,2),'NumLags',winCorr*Fs);
            else
                % Auto-correlation
                [autocorr_tmp1,lags_tmp1]=autocorr(fp_sub_samehemi(:,1),'NumLags',winCorr*Fs);
                [autocorr_tmp2,lags_tmp2]=autocorr(fp_sub_samehemi(:,2),'NumLags',winCorr*Fs);
            end
    
            % Autocorr for DA-ACh pairs channels (saves ACh from autocorr_tmp1 and DA from autocorr_tmp2)
            autocorr_cell{1,1}(:,1)=lags_tmp1;
            autocorr_cell{1,1}(:,2*counter)=autocorr_tmp1;
            autocorr_cell{1,1}(:,2*counter+1)=autocorr_tmp2;
        end

    counter=counter+1;
    end

    waitbar(x/length(beh),h);
end
close(h);

%% AVERAGE ACROSS ALL RECORDINGS FOR ONE ANIMAL SUCH THAT N = X mice


    % Corr across hemisphere (i.e. DA vs DA and ACh vs ACh)
    rec = {}; 
    for x = 1:length(beh)*2; 
        if (mod(x,2)~=0) 
            rec{x} = strcat(strtok(beh(ceil(x/2)).rec,'-'),'ACh'); 
        else 
            rec{x} = strcat(strtok(beh(ceil(x/2)).rec,'-'),'DA'); 
        end 
    end
    uni = unique(rec); 
    
    for i=1:length(rec)
        rec_miceID{i}=rec{i}(1:5);
    end
    
    for i=1:length(uni)
        uni_miceID{i}=uni{i}(1:5);
    end
    
    ACh_rec=contains(rec,"ACh");
    DA_rec=contains(rec,"DA");
    
    uni_miceID=unique(uni_miceID);
    nAn = length(uni_miceID); % number of unique animal IDs (Not the number of unique channels include per animal like DAvsDA)

% For Normal Xcorr
if (channel_type == 2)
    rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end
    uni = unique(rec); 
    nAn = length(uni); % number of unique animal IDs
end

corr_achda = cell(3,1); % initialize output
corr_AChACh = cell(3,1); % initialize output
corr_DADA = cell(3,1);  % initialize output
shuff_achda = cell(3,3); % initialize output 
shuffAChACh = cell(3,3); % initialize output
shuffDADA = cell(3,3);  % initialize output

% For normal
% for x = 1:nAn
%     idx = strcmp(rec,uni{x}); % match animal ID to recordings
%     for z = 1:nStates % iterate over behavioral states
%         corr_adj = corr_cell{z,1}; % extract cross-correlation output for this behavioral state
%         corr_adj = corr_adj - nanmean(corr_adj(1:find(lags./Fs == -2),:)); % adjust such that baseline outside of +/- 2s is at zero
%         corr_achda{z,1}(:,x) = nanmean(corr_adj(:,idx),2); % average across all recordings for this animal
%         for b = 2:4 % iterate over shuffle percentiles
%             corr_adj = corr_cell{z,b};
%             shuff_achda{z,b-1} = nanmean(corr_adj(:,idx),2); % average across all recordings for this animal
%         end
%     end
% end

% Initialize autocorrelation cell array for each channel type
if (behavior_type==1)
    autocorr_ACh = cell(3,1); % initialize output
    autocorr_DA = cell(3,1); % initialize output 
    for z = 1:nStates
       autocorr_ACh{z,1}=autocorr_cell{z,1}(:,1)/Fs;
       autocorr_DA{z,1}=autocorr_cell{z,1}(:,1)/Fs;
    end
else
    autocorr_ACh = cell(1,1); % initialize output
    autocorr_DA = cell(1,1); % initialize output 
end

% Corr across hemisphere (i.e. DA vs DA and ACh vs ACh)
for x = 1:nAn
    idx = strcmp(rec_miceID,uni_miceID{x}); % match channel ID (within each animal) to recordings
    if (channel_type == 2)
    idx2 = strcmp(rec,uni{x}); % match animal ID to recordings. Each index corresponds to animal.
    end
      if (behavior_type==1)
        for z = 1:nStates % iterate over behavioral states
            % Indexing: Gets the corresponding DA or ACh 
            [cACh,~,~]=intersect(find(double(idx)),find(double(ACh_rec)));
            [cDA,~,~]=intersect(find(idx),find(DA_rec));
    
            if (analysis_type==2)
                corr_adj_acrosshemi = corr_cell_acrosshemi{z,1}; % extract cross-correlation output for this behavioral state
                % Idea here is that correlation across different channels/recordings are averaged relative to the baseline at -2
                % seconds in this case, across all corresponding rec#_channelcomparision#
                % For example, above is done across all DAACH xcorr values for all recordings.
                % For DAvsDA and AChvsACh, these 2 xcorr channels must be separated across all recordings.
                % Adjust such that baseline outside of +/- 2s is at zero (for every recording- i.e. DAACh)
                %corr_adj = corr_adj - nanmean(corr_adj(1:find(lags./Fs == -2),:)); 
                % Adjust such that the baseline outside of +/2 is at zero (ONLY for AChACh)
                %corr_adjAChACh = corr_adj - nanmean(corr_adj(1:find(lags./Fs == -2),mod(1:size(corr_adj,2),2)~=0)); % hardcoded
                corr_adjAChACh = corr_adj_acrosshemi - nanmean(corr_adj_acrosshemi(1:find(lags./Fs == -2),ACh_rec)); % non-hardcoded
                % Adjust such that the baseline outside of +/2 is at zero (ONLY for DADA)
                %corr_adjDADA = corr_adj - nanmean(corr_adj(1:find(lags./Fs == -2),mod(1:size(corr_adj,2),2)==0)); %hardcoded
                corr_adjDADA = corr_adj_acrosshemi - nanmean(corr_adj_acrosshemi(1:find(lags./Fs == -2),DA_rec)); % non-hardcoded
                %corr_achda{z,1}(:,x) = nanmean(corr_adj(:,idx),2); % average across all recordings for this animal
                corr_AChACh{z,1}(:,x) = nanmean(corr_adjAChACh(:,cACh),2); % average across all recordings for this animal
                corr_DADA{z,1}(:,x) = nanmean(corr_adjDADA(:,cDA),2); % average across all recordings for this animal
            end
    
                corr_adj_samehemi = corr_cell_samehemi{z,1}; % extract cross-correlation output for this behavioral state
                % Idea here is that correlation across different channels/recordings are averaged relative to the baseline at -2
                % seconds in this case, across all corresponding rec#_channelcomparision#
                % For example, above is done across all DAACH xcorr values for all recordings.
                % For DAvsDA and AChvsACh, these 2 xcorr channels must be separated across all recordings.
                % Adjust such that baseline outside of +/- 2s is at zero (for every recording- i.e. DAACh)
                %corr_adj = corr_adj - nanmean(corr_adj(1:find(lags./Fs == -2),:)); 
                corr_adj_samehemi = corr_adj_samehemi - nanmean(corr_adj_samehemi(1:find(lags./Fs == -2),:)); % adjust such that baseline outside of +/- 2s is at zero
            if (analysis_type == 2)   
                corr_achda{z,1}(:,x) = nanmean(corr_adj_samehemi(:,idx),2); % average across all recordings for this animal
            else
                corr_achda{z,1}(:,x) = nanmean(corr_adj_samehemi(:,idx2),2); % average across all recordings for this animal
            end
    
    
    
           % Extract autocorrelation output for this behavioral state and separate them into eawch channel type (i.e. DA and ACh):
           aACh=[];
           aDA=[];
           if (channel_type == 1) 
              for i=1:length(cACh)
                aACh=[aACh cACh(i)*2 cACh(i)*2+1];
                aDA=[aDA cDA(i)*2 cDA(i)*2+1];
              end
           else
              for i=1:length(cACh)
                aACh=[aACh cACh(i)+1];
                aDA=[aDA cDA(i)+1];
              end
           end
            autocorr_adjACh = autocorr_cell{z,1}(:,aACh);
            autocorr_adjDA = autocorr_cell{z,1}(:,aDA);
            % Save the output
            autocorr_ACh{z,1}(:,x+1) = nanmean(autocorr_adjACh,2);
            autocorr_DA{z,1}(:,x+1) = nanmean(autocorr_adjDA,2);
    
    
            for b = 2:4 % iterate over shuffle percentiles
                if (analysis_type == 2)
                    corr_adj_acrosshemi = corr_cell_acrosshemi{z,b};
                    % Assuming, we don't need to average to the baseline of -2s
                    % because this is a shuffled percentile
                    %[c,ia,ib]=intersect(find(idx),find(ACh_rec))
                    shuffAChACh{z,b-1}(:,x) = nanmean(corr_adj_acrosshemi(:,cACh),2);
                    shuffDADA{z,b-1}(:,x) = nanmean(corr_adj_acrosshemi(:,cDA),2);
                    %shuffAChACh{z,b-1}(:,x) = nanmean(corr_adj(:,ismember(find(idx),find(ACh_rec))),2);
                    %shuffDADA{z,b-1}(:,x) = nanmean(corr_adj(:,ismember(find(idx),find(DA_rec))),2);
                    %shuff_achda{z,b-1} = nanmean(corr_adj(:,idx),2); % average across all recordings for this animal
                end
                corr_adj_samehemi = corr_cell_samehemi{z,b};
                if (channel_type == 1)
                    shuff_achda{z,b-1} = nanmean(corr_adj_samehemi(:,idx),2); % average across all recordings for this channel ID within animal
                else
                    shuff_achda{z,b-1} = nanmean(corr_adj_samehemi(:,idx2),2); % average across all recordings for this animal
                end
            end 
        end
      else
        % Indexing: Gets the corresponding DA or ACh 
        [cACh,~,~]=intersect(find(double(idx)),find(double(ACh_rec)));
        [cDA,~,~]=intersect(find(idx),find(DA_rec));

        if (analysis_type==2)
            corr_adj_acrosshemi = corr_cell_acrosshemi{1,1}; % extract cross-correlation output for this behavioral state
            % Idea here is that correlation across different channels/recordings are averaged relative to the baseline at -2
            % seconds in this case, across all corresponding rec#_channelcomparision#
            % For example, above is done across all DAACH xcorr values for all recordings.
            % For DAvsDA and AChvsACh, these 2 xcorr channels must be separated across all recordings.
            % Adjust such that baseline outside of +/- 2s is at zero (for every recording- i.e. DAACh)
            %corr_adj = corr_adj - nanmean(corr_adj(1:find(lags./Fs == -2),:)); 
            % Adjust such that the baseline outside of +/2 is at zero (ONLY for AChACh)
            %corr_adjAChACh = corr_adj - nanmean(corr_adj(1:find(lags./Fs == -2),mod(1:size(corr_adj,2),2)~=0)); % hardcoded
            corr_adjAChACh = corr_adj_acrosshemi - nanmean(corr_adj_acrosshemi(1:find(lags./Fs == -2),ACh_rec)); % non-hardcoded
            % Adjust such that the baseline outside of +/2 is at zero (ONLY for DADA)
            %corr_adjDADA = corr_adj - nanmean(corr_adj(1:find(lags./Fs == -2),mod(1:size(corr_adj,2),2)==0)); %hardcoded
            corr_adjDADA = corr_adj_acrosshemi - nanmean(corr_adj_acrosshemi(1:find(lags./Fs == -2),DA_rec)); % non-hardcoded
            %corr_achda{z,1}(:,x) = nanmean(corr_adj(:,idx),2); % average across all recordings for this animal
            corr_AChACh{1,1}(:,x) = nanmean(corr_adjAChACh(:,cACh),2); % average across all recordings for this animal
            corr_DADA{1,1}(:,x) = nanmean(corr_adjDADA(:,cDA),2); % average across all recordings for this animal
        end

            corr_adj_samehemi = corr_cell_samehemi{1,1}; % extract cross-correlation output for this behavioral state
            % Idea here is that correlation across different channels/recordings are averaged relative to the baseline at -2
            % seconds in this case, across all corresponding rec#_channelcomparision#
            % For example, above is done across all DAACH xcorr values for all recordings.
            % For DAvsDA and AChvsACh, these 2 xcorr channels must be separated across all recordings.
            % Adjust such that baseline outside of +/- 2s is at zero (for every recording- i.e. DAACh)
            %corr_adj = corr_adj - nanmean(corr_adj(1:find(lags./Fs == -2),:)); 
            corr_adj_samehemi = corr_adj_samehemi - nanmean(corr_adj_samehemi(1:find(lags./Fs == -2),:)); % adjust such that baseline outside of +/- 2s is at zero
        if (analysis_type == 2)   
            corr_achda{1,1}(:,x) = nanmean(corr_adj_samehemi(:,idx),2); % average across all recordings for this animal
        else
            corr_achda{1,1}(:,x) = nanmean(corr_adj_samehemi(:,idx2),2); % average across all recordings for this animal
        end
    
       % Extract autocorrelation output for this behavioral state and separate them into eawch channel type (i.e. DA and ACh):
       aACh=[];
       aDA=[];
       if (channel_type == 1) 
          for i=1:length(cACh)
            aACh=[aACh cACh(i)*2 cACh(i)*2+1];
            aDA=[aDA cDA(i)*2 cDA(i)*2+1];
          end
       else
          for i=1:length(cACh)
            aACh=[aACh cACh(i)+1];
            aDA=[aDA cDA(i)+1];
          end
       end
        autocorr_adjACh = autocorr_cell{1,1}(:,aACh);
        autocorr_adjDA = autocorr_cell{1,1}(:,aDA);
        % Save the output
        autocorr_ACh{1,1}(:,x+1) = nanmean(autocorr_adjACh,2);
        autocorr_DA{1,1}(:,x+1) = nanmean(autocorr_adjDA,2);


        for b = 2:4 % iterate over shuffle percentiles
            if (analysis_type == 2)
                corr_adj_acrosshemi = corr_cell_acrosshemi{1,b};
                % Assuming, we don't need to average to the baseline of -2s
                % because this is a shuffled percentile
                %[c,ia,ib]=intersect(find(idx),find(ACh_rec))
                shuffAChACh{1,b-1}(:,x) = nanmean(corr_adj_acrosshemi(:,cACh),2);
                shuffDADA{1,b-1}(:,x) = nanmean(corr_adj_acrosshemi(:,cDA),2);
                %shuffAChACh{z,b-1}(:,x) = nanmean(corr_adj(:,ismember(find(idx),find(ACh_rec))),2);
                %shuffDADA{z,b-1}(:,x) = nanmean(corr_adj(:,ismember(find(idx),find(DA_rec))),2);
                %shuff_achda{z,b-1} = nanmean(corr_adj(:,idx),2); % average across all recordings for this animal
            end
            corr_adj_samehemi = corr_cell_samehemi{1,b};
            if (channel_type == 1)
                shuff_achda{1,b-1} = nanmean(corr_adj_samehemi(:,idx),2); % average across all recordings for this channel ID within animal
            else
                shuff_achda{1,b-1} = nanmean(corr_adj_samehemi(:,idx2),2); % average across all recordings for this animal
            end
        end
      end


end

end