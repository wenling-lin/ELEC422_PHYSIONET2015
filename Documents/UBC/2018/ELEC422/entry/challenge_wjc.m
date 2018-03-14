%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ELEC 422 - Physionet 2015 Challenge
% Authors:  Wen-Ling Lin (31120132)
%           Janelle Somerville (52155132)
%           Claire Chang (35340132)
%
% Sources Used: 
% Pre-Processing Filter Specs Used: http://www.cinc.org/archives/2015/pdf/1181.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function alarmResult=challenge(recordName,alarm_type)
function alarmResult =challenge_(recordName, alarm_type);
%clear; clear all; clc;
answers = 'answers.txt';
% recordName = 'a311l'; Good Data across all channels

% recordName = 'b285l'; Good PPG 
% recordName = 'b286s'; %-- OMG BREAKS IT
% recordName = 'b220s'; % Has PPG Signal but clipped -> ECG reliable
% recordName = 'b216s'; % Has clean PPG signal
% recordName = 'a442s'; % True Asystole Alarm
% recordName = 'a429l'; % False Asystole Alarm
% recordName = 'a103l'; % False Alarm
% recordName = 'a302s'; % False Alarm
% recordName = 'a673l'; % False Alarm
% recordName = 'a639l'; % True Alarm
% recordName = 'a186s';
%recordName = 'a272s'; % False Alarm
%recordName = 'a152s'; % Why this no work.... Whats wrong with wavedec
% recordName = 'a167l';
% recordName = 'a171l';
% recordName = 'a391l'; % One peak but top of a curve (ECG) false asystole
% recordName = 'a443l'; % One peak but true asystole

%%
%Get all ECG, blood pressure and photoplethysmogram signals
[~,signal,Fs,siginfo]=rdmat(recordName);
alarmResult=1;
description=squeeze(struct2cell(siginfo));
description=description(4,:);

% Resample signal to 125Hz
Fs=Fs(1);
if Fs~=125
    signal=resample(signal,125,Fs);
    Fs=125;
end

%% Aquire Signals and Signal Quality Index - Ensure Signals are Viable Before Processing 
%   Run WABP on the record, which by default will analyze the first ABP, ART, or BP signal
% Make sure signal does not have any "NaN"s
% ECG II
ecg1_ind = get_index(description, 'II');
ecg1 = signal(:,ecg1_ind);
ecg1 = ecg1(~isnan(ecg1));

% ECG V
ecg2_ind = get_index(description,'V');
ecg2 = signal(:,ecg2_ind);

ecg_signal = ecg1';

% set valid data segment for decision making, 16s before the alarm
%N_d=Fs*5*60; % alarm position % 5 minute point
%N0_d=N_d-Fs*16+1; % 16s before the alarm % 16 seconds before 5 minute point

N_d=Fs*5*60; % alarm position % 5 minute point
N0_d=N_d-Fs*16; % 16s before the alarm % 16 seconds before 5 minute point

if(length(ecg1) < N_d)
    alarmResult = 'False Alarm'
    alarmResult = 0;
    return;
end

%ecg_signal = ecg_signal(ecg_signal != NaN);
[C,L]=wavedec(ecg_signal,8,'db4'); % [approx, detail] as per lecture
[thr,sorh,keepapp]=ddencmp('den','wv',ecg_signal);
cleanecg=wdencmp('gbl',C,L,'db4',8,thr,sorh,keepapp);

% Clean unneeded variables - Prevents the MEMORY Error from occuring
clear f_y; clear C; clear L;

% Remove Baseline Wander
opol = 6; %polynomial degree
ecg_signal = ecg1';
[p,s,mu] = polyfit((1:numel(ecg_signal)),ecg_signal,opol);
f_y = polyval(p,(1:numel(ecg_signal)),[],mu);
ecg_signal_dt = ecg_signal - f_y;


N=[];
N0=[];
abp_ind=get_index(description,'ABP');
ann_abp=[];
features=[];
BEATQ=[];
R=[];
if(~isempty(abp_ind))
   % Pre-Filter ABP Data using Bandpass with Cut-off = 0.5 and 10 Hz
   abp = signal(:,abp_ind);
   [b,a] = butter(3,[0.5 10]/(Fs/2),'bandpass');
   abp_f = filter(b,a,abp);
   ann_abp=wabp(abp,0,1);
   
   % Analyze the signal quality index of ABP using jSQI
   if length(ann_abp)>=3 % at least 3 abp beats detected
        [features] = abpfeature(signal(:,abp_ind),ann_abp);
        [BEATQ R] = jSQI(features, ann_abp, signal(:,abp_ind));
   end
end

%Run WABP on the record of 'PLETH' to analyze photoplethysmogram signal
ppg_ind=get_index(description,'PLETH');
ann_ppg=[];
if (~isempty(ppg_ind))
    % Pre-Filter PLETH Data using Bandpass with Cut-off = 0.5 and 10 Hz
    ppg = signal(:, ppg_ind);
    [b,a] = butter(3,[0.5 10]/(Fs/2),'bandpass');
    ppg_f = filter(b,a,ppg);
   
    y=quantile(signal(:,ppg_ind),[0.05,0.5,0.95]);
    ann_ppg=wabp(signal(:,ppg_ind),0,(y(3)-y(1))/120);
    % Analyze the signal quality index of PPG 
    if ~isempty(ann_ppg)
        [psqi]=ppgSQI(signal(:,ppg_ind),ann_ppg);
    end
end

% % set valid data segment for decision making, 16s before the alarm
% N_d=Fs*5*60; % alarm position % 5 minute point
% N0_d=N_d-Fs*16+1; % 16s before the alarm % 16 seconds before 5 minute point

% if(length(ecg1) < N_d)
%     alarmResult = 'False Alarm'
%     alarmResult = 0;
%     return;
% end
%% Find the Beats and calculate Heart Rate from each Channel
% select the beats in the segment
n_abp_beats=intersect(find(ann_abp>=N0_d),find(ann_abp<=N_d));
n_ppg_beats=intersect(find(ann_ppg>=N0_d),find(ann_ppg<=N_d));

hr_max_abp=NaN;
max_rr_abp=NaN;
hr_max_ppg=NaN;
max_rr_ppg=NaN;
ecg_minpkdistance=NaN;

pk2pk_ppg_HR = NaN;
pk2pk_abp_HR = NaN;
pk2pk_ecg_HR = NaN;

all_channel_HR = [];

% ABP
% calculate the heart rate using ABP - wabp beat to beat
if length(n_abp_beats)>=2
    hr_max_abp=60*Fs/min(diff(ann_abp(n_abp_beats)));
    max_rr_abp=max(diff(ann_abp(n_abp_beats)))/Fs;
    
    % Try peak detection method using the filtered data
    % Restrict to window of 16 seconds before alarm
    abp_f = abp_f(N0_d:N_d);
    [abp_pks,abp_locs] = findpeaks(abp_f, 'MinPeakHeight', 1.25* mean(abs(abp_f)));
    findpeaks(abp_f, 'MinPeakHeight', 1.25 * mean(abs(abp_f)));
    
    if( strcmp(alarm_type,'Asystole') )
        abp_locs = vertcat(abp_locs, (N_d - N0_d + 1));
    end
    
    ecg_minpkdistance = min(abp_locs); % Use this minimum peak distance later for ecg peak detector    
    hr_max_abp_pk2pk = 60*Fs/min(diff(abp_locs));
    pk2pk_abp_HR = 60*Fs ./ diff(abp_locs)
    % Determine standard deviation and remove heart rate +/- less than 
    % 3 standard deviations from the mean
    abp_std = std(pk2pk_abp_HR);
    abp_median = median(pk2pk_abp_HR);
    pk2pk_abp_HR = pk2pk_abp_HR((pk2pk_abp_HR > 0) & (pk2pk_abp_HR < (abp_median +(3 * abp_std))))
    
    % Calculate the HR
    % Physionet
    abp_phys_hr=60*Fs ./diff(ann_abp(n_abp_beats));
    abp_phys_mean = mean(abp_phys_hr)
    % Pk2PK
    abp_pk2pk_mean = mean(pk2pk_abp_HR)
    % Add all HR into the channel median vector
    all_channel_HR = vertcat(all_channel_HR,abp_phys_hr, pk2pk_abp_HR);
end

% PPG
% calculate the heart rate using PPG 
if length(n_ppg_beats)>=2
    hr_max_ppg=60*Fs/min(diff(ann_ppg(n_ppg_beats)));
    max_rr_ppg=max(diff(ann_ppg(n_ppg_beats)))/Fs;
    ppg_HR = 60*Fs ./ diff(ann_ppg(n_ppg_beats));
    
    % Try peak detection method using the filtered data
    % Restrict to window of 16 seconds before alarm
    ppg_f = ppg_f(N0_d:N_d);
    [ppg_pks, ppg_locs] = findpeaks(ppg_f, 'MinPeakHeight', 1.25 * mean(abs(ppg_f)));
    findpeaks(ppg_f, 'MinPeakHeight', 1.25 * mean(abs(ppg_f)))
    ecg_minpkdistance = min(ppg_locs); % Use this minimum peak distance later for ecg peak detector
    
    if( strcmp(alarm_type,'Asystole') )
        ppg_locs = vertcat(ppg_locs, (N_d - N0_d + 1));
    end
    
    hr_max_ppg_pk2pk = 60*Fs/min(diff(ppg_locs));
    pk2pk_ppg_HR = 60*Fs ./ diff(ppg_locs);
    % Determine standard deviation and remove heart rate +/- less than 
    % 3 standard deviations from the mean
    ppg_std = std(pk2pk_ppg_HR);
    ppg_median = median(pk2pk_ppg_HR);
    pk2pk_ppg_HR = pk2pk_ppg_HR((pk2pk_ppg_HR > 0) & (pk2pk_ppg_HR < (ppg_median +(3 * ppg_std)))) 

    % Physionet
    ppg_phys_hr=60*Fs ./diff(ann_ppg(n_ppg_beats));
    ppg_phys_mean = mean(ppg_phys_hr)
    % Pk2PK
    ppg_pk2pk_mean = mean(pk2pk_ppg_HR)
    % Add all HR into the channel median vector
    all_channel_HR = vertcat(all_channel_HR,ppg_phys_hr, pk2pk_ppg_HR)
end

% ECG
% calculate the heart rate using ECG
% Restrict to window of 16 seconds before alarm
    % Assess each channel the derived heart rate to choose
    % between minpeakheight or minpeakdistance for pk detection
    % Define 30-240 as beats !! Would probably need to modify for Asystole
   
    ecg1 = ecg_signal_dt(N0_d:N_d);
    % Flag to say that the HR has been calc
    hr_ecg_calculated = 0;
    if((~isempty(abp_ind) & length(n_abp_beats)>=2 | isnan(ecg_minpkdistance)))
        if(length(pk2pk_abp_HR < 30 | pk2pk_abp_HR > 240) > round((length(pk2pk_abp_HR)/4), 0) )
            % Enter this condition = ABP channel clean
            [ecg1_pks,ecg1_locs] = findpeaks(ecg1, 'MinPeakHeight', 1.25*max(ecg1)/2);
            findpeaks(ecg1, 'MinPeakHeight', 1.25*max(ecg1)/2);
            hr_ecg_calculated = 1;
        end
    end
    if(~isempty(ppg_ind) & length(n_ppg_beats)>=2 | isnan(ecg_minpkdistance) )
        if(length(pk2pk_ppg_HR < 30 | pk2pk_ppg_HR > 240) > round((length(pk2pk_ppg_HR)/4), 0) )
            % Enter this condition = PPG channel clean
            [ecg1_pks,ecg1_locs] = findpeaks(ecg1, 'MinPeakHeight', max(ecg1)/2);
            findpeaks(ecg1, 'MinPeakHeight', max(ecg1)/2);
            hr_ecg_calculated = 1;
        end
    end
    % This will calculate using the minPeakDistance by calculating the min
    % samples from the other channels. If the data is noisy from PPG or ABP
    if( ~hr_ecg_calculated )
        [ecg1_pks,ecg1_locs] = findpeaks(ecg1, 'MinPeakDistance', 1.25*ecg_minpkdistance );
        findpeaks(ecg1, 'MinPeakDistance', 1.25*ecg_minpkdistance );
        %findpeaks(ecg1, 'MinPeakHeight', max(ecg1)/2);
    end

    if( length(ecg1_locs) > 1 )
        
        if( strcmp(alarm_type, 'Asystole') )
            ecg1_locs = horzcat(ecg1_locs, (N_d - N0_d + 1));
        end
        
        hr_max_ecg_pk2pk = 60*Fs/min(diff(ecg1_locs));
        max_rr_ecg=max(diff(ecg1_locs))/Fs;
        pk2pk_ecg_HR = 60*Fs ./ diff(ecg1_locs);
        pk2pk_ecg_HR = pk2pk_ecg_HR';
        % Determine standard deviation and remove heart rate +/- less than 
        % 3 standard deviations from the median (median robust to extreme HR)
        ecg_std = std(pk2pk_ecg_HR);
        ecg_median = median(pk2pk_ecg_HR);
        pk2pk_ecg_HR = pk2pk_ecg_HR((pk2pk_ecg_HR > 0) & (pk2pk_ecg_HR < (ecg_median +(3 * ecg_std)))) 
        % Calculate HR
        ecg_pk2pk_mean = mean(pk2pk_ecg_HR)
        % Add all HR into the channel median vector
        all_channel_HR = vertcat(all_channel_HR,pk2pk_ecg_HR)
    end
    
%% calculate the signal quality index
if ~isempty(ann_abp)
    abpsqi=1-sum(sum(BEATQ(intersect(n_abp_beats,1:length(BEATQ)),:)))/numel(BEATQ(intersect(n_abp_beats,1:length(BEATQ)),:));
else
    abpsqi=0;
end
if ~isempty(ann_ppg)
    ppgsqi=mean(psqi(intersect(n_ppg_beats,1:length(psqi))));
else
    ppgsqi=0;
end

% Signal Quality Index threshold
sqi_th = 0.9;
abpsqi;
ppgsqi;

%% Decide which channel has the cleanest signal to select the HR from 
% 1. Calculate the HR in all channels and concatenate into one big vector
% 2. Calculate the median of the big vector of combined HRs from all channels
% 3. Select Channel that has the closest median to the HR from all channels

% Calculate median of all calculated heart rates from all channels
% Too sensitive to high outliers so use median
all_channel_median = median(all_channel_HR<300)
selected_signal = [];
% Compare means (Using mean because want as few outliers from the 
% signal in the selected channel

% Temporary Variables - 0 or Nan (Decide Later)
smallest_delta = 999999999; %Just arbitrarily set a large delta
selected_channel = [];
% ABP 
if (~isempty(abp_ind) & length(n_abp_beats)>=2)
    physio_delta = abs(abp_phys_mean - all_channel_median);
    pk2pk_delta = abs(abp_pk2pk_mean - all_channel_median);
    
    % Physionet is closer to mean
    if (physio_delta < pk2pk_delta )
        % Compare with whats been caculated (Physio is smallest) ABP
        if(physio_delta < smallest_delta & abpsqi >= sqi_th) % Added the & !!! - March 12
            selected_signal = abp_phys_hr;
            smallest_delta = physio_delta;
            selected_channel = 'ABP';
        end
    else
        % Pk2pk is smallest ABP
        if(pk2pk_delta < smallest_delta)
            selected_signal = pk2pk_abp_HR;
            smallest_delta = pk2pk_delta;
            selected_channel = 'ABP';
        end
    end
end

% PPG
if(~isempty(ppg_ind) & length(n_ppg_beats)>=2 )
    physio_delta = abs(ppg_phys_mean - all_channel_median);
    pk2pk_delta = abs(ppg_pk2pk_mean - all_channel_median);
    
        % Physionet is closer to mean
    if (physio_delta < pk2pk_delta )
        % Compare with whats been caculated (Physio is smallest) PPG
        if(physio_delta < smallest_delta)
            selected_signal = ppg_phys_hr;
            smallest_delta = physio_delta;
            selected_channel = 'PPG';
        end
    else
        % Pk2pk is smallest PPG
        if(pk2pk_delta < smallest_delta & ppgsqi >= sqi_th) % Added the %
            selected_signal = pk2pk_ppg_HR;
            smallest_delta = pk2pk_delta;
            selected_channel = 'PPG';
        end
    end
end

% ECG
if(length(ecg1_locs) > 1)
    pk2pk_delta = abs(ecg_pk2pk_mean - all_channel_median);
    
    if(pk2pk_delta < smallest_delta)
        selected_signal = pk2pk_ecg_HR;
        smallest_delta = pk2pk_delta;
        selected_channel = 'ECG';
    end
end

% Extra Signal Check -- If ECG Channel is empty then take from abp or ppg 
if isempty(selected_channel)
    if( abpsqi ~= 0 )
        selected_channel = 'abp';
        selected_signal = pk2pk_abp_HR;
    elseif( ppgsqi ~=0)
        selected_channel = 'ppg';
        selected_signal = pk2pk_ppg_HR;
    end
end
%%  Signal from cleanest channel has been selected - Calculate Heart Rate and Detect Arrythmia Type

smallest_delta
selected_signal
% Recalculate hr_max and max_rr
% Selected signal will be empty if in the window of samples there is only
% one peak. In the case of asystole this could mean that the person has
% flatlined. Note: This is an edge case but should account for


% If only one peak and asystole alarm
if( (length(ecg1_locs) == 1) & (alarm_type == 'Asystole'))
        alarmResult = 'True Alarm'
        alarmResult = 1;
        return;
end
        
hr_max=60*Fs/min(selected_signal);
max_rr=max(selected_signal)/Fs;

%% Bradycardia
% Redefine -> Under 50 beats
% 3 Consecutive Less than 50
brady_thresh = 50;
isBrady = 0;
brady_beats = 3;
if length(selected_signal) >= brady_beats
    brady_count = 0;
    for i = 1:length(selected_signal)
        if( selected_signal(i) < brady_thresh )
            brady_count = brady_count + 1;
            if brady_count == brady_beats
                isBrady = 1;
                break;
            end
        else
            brady_count = 0;
        end
    end
end

%% Decision Making -- Suppress or Don't Suppress Alarm 

% Alarm threshold (seconds)
ASY_th = 4;
BRA_th = 40;
TAC_th = 140;
VTA_th = 100;
VFB_th = 250;
tolerance = 5; % tolerance = 5 bmp

% if(strncmpi(recordName,'a',1))
%     alarm_type = 'a';
% end
switch alarm_type
    case 'Asystole'
    %case 'a'
        % No QRS for at least 4 seconds
        if (max_rr<ASY_th && sum(selected_signal < 15) == 0)
            alarmResult = 'False Alarm'
            alarmResult = 0;
        else
            alarmResult = 'True Alarm'
            alarmResult = 1;
        end
    case 'Bradycardia'
        if (isBrady == 1)            
            alarmResult = 'True Alarm'
            alarmResult = 1;
        else
            alarmResult = 'False Alarm'
            alarmResult = 0;
        end
    otherwise
        error(['Unknown alarm type: ' alarm_type])
end

% Write result to answers.txt
fid = fopen(answers, 'a');
if (fid == -1)
    error('Could not open answer file');
end

% Get base name of record (without directories)
i = strfind(recordName, filesep);
if (~isempty(i))
    basename = recordName(i(end)+1 : end);
else
    basename = recordName;
end

fprintf(fid, '%s,%d\n', basename, alarmResult);
fclose(fid);
end
%% Plotting

%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%
function ind=get_index(description,pattern)
M=length(description);
tmp_ind=strfind(description,pattern);
ind=[];
for m=1:M
    if(~isempty(tmp_ind{m}))
        ind(end+1)=m;
    end
end
end
