function [corr_achda, lags, shuff_achda] = AK_corrFP(beh, varargin)
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

%% OUTPUTS
% mat = struct; % temporary output structure
corr_cell = cell(3,4); % temporary output cell array
for a = 1:3; for b = 1:4; corr_cell{a,b} = nan((winCorr*2*50)+1,length(beh)); end; end % fill cell array with nan's


%% RUN ANALYSIS ON ALL RECORDINGS
h = waitbar(0, 'cross correlation');
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

    % Corr across hemisphere (i.e. DA vs DA and ACh vs ACh)
   for xx = 1:length(beh(x).FP)/2
     y = [xx+2 xx];
     fp_mat = [];
     fp_mat(:,1) = beh(x).FP{y(1)}; Fs = beh(x).Fs; % extract photometry signal from structure, which will be used as reference
     fp_mat(:,1) = fp_mat(:,1) - nanmean(fp_mat(:,1)); % subtract baseline (mean of entire photometry signal) from fp
     fp_mat(:,2) = beh(x).FP{y(2)}; % extract photometry signal from structure
     fp_mat(:,2) = fp_mat(:,2) - nanmean(fp_mat(:,2)); % subtract baseline (mean of entire photometry signal) from fp



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
    fp_sub = [];
    for z = 1:nStates
        if length(idxStates{x,z})< 2; continue; end % Only calculate states that have more than 2 indicies
        fp_sub = fp_mat(idxStates{x,z},:); % signal
        
        % fp_sub = normalize(fp_mat,1,'range'); % normalize [0 1]
        % fp_sub = fp_filt(idx_cell{z},:); % bandpass filtered signal
        % [xcf, lags, bounds] = crosscorr(fp_sub(:,1), fp_sub(:,2),'NumLags',100,'NumSTD',3);
        % [shuff,~,~] = crosscorr(fp_sub(randperm(size(fp_sub,1)),1), fp_sub(randperm(size(fp_sub,2)),2),'NumLags',100,'NumSTD',3);

        % Auto-correlation
        [autocorr_tmp1,lags_tmp1]=autocorr(fp_sub(:,1),'NumLags',winCorr*Fs);
        [autocorr_tmp2,lags_tmp2]=autocorr(fp_sub(:,2),'NumLags',winCorr*Fs);

        % Cross-Correlation
        [corr_tmp, lags] = xcorr(fp_sub(:,1), fp_sub(:,2), winCorr*Fs, 'coeff'); % cross-correlation
        fp_sub_new = fp_sub(:,2);
        tmp_shuff = []; 
        for s = 1:50
            fp_sub_new = circshift(fp_sub_new, Fs);
            % tmp_shuff(:,s) = xcorr(fp_sub(randperm(size(fp_sub,1)),1), fp_sub(randperm(size(fp_sub,2)),2), 10*Fs, 'coeff');
            % tmp_shuff(:,s) = xcorr(fp_sub(:,1), fp_sub(randperm(size(fp_sub,2)),2), 10*Fs, 'coeff');
            tmp_shuff(:,s) = xcorr(fp_sub(:,1), fp_sub_new, winCorr*Fs, 'coeff');
        end

        % Normal DA vs ACh 2-channel xcorr
        %corr_cell{z,1}(:,x) = corr_tmp;                     % cross-correlation
        %corr_cell{z,2}(:,x) = prctile(tmp_shuff, 5, 2);     % shuffle 5th percentile
        %corr_cell{z,3}(:,x) = prctile(tmp_shuff, 50, 2);    % shuffle 50th percentile
        %corr_cell{z,4}(:,x) = prctile(tmp_shuff, 95, 2);    % shuffle 95th percentile

        % Corr across hemisphere (i.e. DA vs DA and ACh vs ACh)
        corr_cell{z,1}(:,counter) = corr_tmp;                  % cross-correlation
        corr_cell{z,2}(:,counter) = prctile(tmp_shuff, 5, 2);  % shuffle 5th percentile
        corr_cell{z,3}(:,counter) = prctile(tmp_shuff, 50, 2); % shuffle 50th percentile
        corr_cell{z,4}(:,counter) = prctile(tmp_shuff, 95, 2); % shuffle 95th percentile
        % Autocorr for 2x2 channels (i.e. DA and DA and ACh and ACh)
        autocorr_cell{z,1}(:,1)=lags_tmp1;
        autocorr_cell{z,1}(:,2*counter)=autocorr_tmp1;
        autocorr_cell{z,1}(:,2*counter+1)=autocorr_tmp2;

    end
    counter=counter+1;
    end

    waitbar(x/length(beh),h);
end
close(h);

%% AVERAGE ACROSS ALL RECORDINGS FOR ONE ANIMAL SUCH THAT N = X mice

% For Normal Xcorr
% rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end
% uni = unique(rec); 
% nAn = length(uni); % number of unique animal IDs

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

corr_achda = cell(3,1); % initialize output
shuff_achda = cell(3,3); % initialize output 

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
autocorr_ACh = cell(3,1); % initialize output
autocorr_DA = cell(3,1); % initialize output 
for z = 1:nStates
   autocorr_ACh{z,1}=autocorr_cell{z,1}(:,1)/Fs;
   autocorr_DA{z,1}=autocorr_cell{z,1}(:,1)/Fs;
end

% Corr across hemisphere (i.e. DA vs DA and ACh vs ACh)
for x = 1:nAn
    idx = strcmp(rec_miceID,uni_miceID{x}); % match animal ID to recordings
    for z = 1:nStates % iterate over behavioral states
        % Indexing: Gets the corresponding DA or ACh 
        [cACh,~,~]=intersect(find(double(idx)),find(double(ACh_rec)));
        [cDA,~,~]=intersect(find(idx),find(DA_rec));

        corr_adj = corr_cell{z,1}; % extract cross-correlation output for this behavioral state
        %autocorr_adj = autocorr_cell_adj{z,1}(:,2:end); % extract the auto-correlation output for this behavioral state
        
        aACh=[];
        aDA=[];
        for i=1:length(cACh)
          aACh=[aACh cACh(i)*2 cACh(i)*2+1];
          aDA=[aDA cDA(i)*2 cDA(i)*2+1];
        end
        autocorr_adjACh = autocorr_cell{z,1}(:,aACh);
        autocorr_adjDA = autocorr_cell{z,1}(:,aDA);


        
        % Idea here is that correlation across different
        % channels/recordings are averaged relative to the baseline at -2
        % seconds in this case, across all corresponding rec#_channelcomparision#
        % For example, above is done across all DAACH xcorr values for all
        % recordings.
        % For DAvsDA and AChvsACh, these 2 xcorr channels must be separated
        % across all recordings.
        % adjust such that baseline outside of +/- 2s is at zero (for every
        % recording- i.e. DAACh)
        %corr_adj = corr_adj - nanmean(corr_adj(1:find(lags./Fs == -2),:)); 
        % Adjust such that the baseline outside of +/2 is at zero (ONLY for AChACh)
        %corr_adjAChACh = corr_adj - nanmean(corr_adj(1:find(lags./Fs == -2),mod(1:size(corr_adj,2),2)~=0)); % hardcoded
        corr_adjAChACh = corr_adj - nanmean(corr_adj(1:find(lags./Fs == -2),ACh_rec)); % non-hardcoded
        
        % Adjust such that the baseline outside of +/2 is at zero (ONLY for DADA)
        %corr_adjDADA = corr_adj - nanmean(corr_adj(1:find(lags./Fs == -2),mod(1:size(corr_adj,2),2)==0)); %hardcoded
        corr_adjDADA = corr_adj - nanmean(corr_adj(1:find(lags./Fs == -2),DA_rec)); % non-hardcoded
        %corr_achda{z,1}(:,x) = nanmean(corr_adj(:,idx),2); % average across all recordings for this animal
        

        corr_AChACh{z,1}(:,x) = nanmean(corr_adjAChACh(:,cACh),2);
        corr_DADA{z,1}(:,x) = nanmean(corr_adjDADA(:,cDA),2);
        %corr_AChACh{z,1}(:,x) = nanmean(corr_adjAChACh(:,ismember(find(idx),find(ACh_rec))),2);
        %corr_DADA{z,1}(:,x) = nanmean(corr_adjDADA(:,ismember(find(idx),find(DA_rec))),2);
        for b = 2:4 % iterate over shuffle percentiles
            corr_adj = corr_cell{z,b};
            % Assuming, we don't need to average to the baseline of -2s
            % because this is a shuffled percentile
            %[c,ia,ib]=intersect(find(idx),find(ACh_rec))
            shuffAChACh{z,b-1}(:,x) = nanmean(corr_adj(:,cACh),2);
            shuffDADA{z,b-1}(:,x) = nanmean(corr_adj(:,cDA),2);
            %shuffAChACh{z,b-1}(:,x) = nanmean(corr_adj(:,ismember(find(idx),find(ACh_rec))),2);
            %shuffDADA{z,b-1}(:,x) = nanmean(corr_adj(:,ismember(find(idx),find(DA_rec))),2);
            %shuff_achda{z,b-1} = nanmean(corr_adj(:,idx),2); % average across all recordings for this animal
        end

      % Extract autocorrelation output for this behavioral state:
       autocorr_ACh{z,1}(:,x+1) = nanmean(autocorr_adjACh,2);
       autocorr_DA{z,1}(:,x+1) = nanmean(autocorr_adjDA,2);
      
    end
end

end