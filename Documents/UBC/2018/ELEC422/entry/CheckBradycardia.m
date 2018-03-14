% Check Bradycardia
clear all;close all;clc
data_dir=[pwd filesep];

%% File Handling

% Write results to the "bradycardia_results.txt" file
fileName = 'bradycardia_results6.txt';
fileID1 = fopen(fileName, 'w');

% Open the alarms true scoring sheet
fid=fopen([data_dir 'ALARMS'],'r');
if(fid ~= -1)
    RECLIST=textscan(fid,'%s %s %d','Delimiter',',');
    fclose(fid);
else
    error('Could not open ALARMS.txt for results. Exiting...')
end

fileName = 'Bradycardia_ALARMS.txt';
fileID2 = fopen(fileName, 'w');

%Add the function on this directory to the MATLAB path
%This is not permanent. This change is only valid for this session of
%MATLAB (will reset once MATLAB is restarted).
addpath(pwd)

if(exist('OCTAVE_VERSION'))
    more off %this seems necessary in order to get back the screen in Octave, but we have not tested this script on Octave yet.
end

%Check for previous files before starting test
answers=dir(['answers.txt']);
fid=fopen([data_dir 'ALARMS'],'r');
if(fid ~= -1)
    RECLIST=textscan(fid,'%s %s %d','Delimiter',',');
    fclose(fid);
else
    error('Could not open ALARMS.txt for scoring. Exiting...')
end

RECORDS=RECLIST{1};
ALARMS=RECLIST{2};
ALARM_RESULTS=RECLIST{3};
N=length(RECORDS);
results=zeros(N,1);
%%
fprintf(fileID1, 'Script ran on: ');
fprintf(fileID1, datestr(now,'HH:MM:SS.FFF\n\n'));
fprintf(fileID2,datestr(now,'HH:MM:SS.FFF\n\n'));

% Only get Bradycardia values:,  
i = contains(RECORDS,'b');
RECORDS = RECORDS(i);
i = contains(ALARMS, 'Bradycardia');
ALARMS = ALARMS(i);
BRADYCARDIA_RESULTS = ALARM_RESULTS(i);

n_correct = 0;
n_falsep = 0;
n_falsen = 0;

for i = 1:length(RECORDS)
    fname=RECORDS{i};
    alarm_scoring = BRADYCARDIA_RESULTS(i);
    
    %results(i)=challenge_wjc([data_dir fname],ALARMS{i});
    results = challenge_wjc([data_dir fname],ALARMS{i})
    
    % Calculate results from our algorithm
    if results == 0
        results = 'False Alarm';
    else
        results = 'True Alarm';
    end
    
    % Parse results from the answer file
    if alarm_scoring == 0
        alarm_scoring = 'False Alarm';
    else
        alarm_scoring = 'True Alarm';
    end
        
    fprintf(fileID1,fname);
    fprintf(fileID1,': ');
    fprintf(fileID1, results);
    fprintf(fileID1, ' => ');
    
    % Filter Scoring Results:
    fprintf(fileID2, fname);
    fprintf(fileID2, ': ');
    fprintf(fileID2, alarm_scoring);
    fprintf(fileID2, '\n');
    
    % Determine the type of True or False 
    if( strcmp(results,alarm_scoring) )
        fprintf(fileID1, 'Correct');
        fprintf(fileID1, '\n');
        n_correct = n_correct + 1;
    elseif( strcmp(results,'True Alarm') & strcmp(alarm_scoring,'False Alarm') )
        fprintf(fileID1, 'False Positive => Calculated True Alarm but Expected False Alarm');
        fprintf(fileID1, '\n');
        n_falsep = n_falsep + 1;
    elseif( strcmp(results, 'False Alarm') & strcmp(alarm_scoring,'True Alarm') )
        fprintf(fileID1, 'False Negative => Calculated False Alarm but Expected True Alarm');
        fprintf(fileID1, '\n');
        n_falsen = n_falsen + 1;
    end
    
    fname
end;

challenge_score = n_correct / (n_correct + n_falsep + 5*n_falsen);

fprintf(fileID1,'=========================================================================\n');
fprintf(fileID1, 'Correct: %d/%d (%f)\n',n_correct, length(RECORDS), 100*(n_correct/length(RECORDS)) );
fprintf(fileID1, 'Incorrect: %d/%d (%f)\n', (n_falsep + n_falsen), length(RECORDS), 100*((n_falsep + n_falsen)/length(RECORDS)) );
fprintf(fileID1, 'False Positive: %d/%d (%f) => Calculated True but Expected False\n', n_falsep, length(RECORDS), 100*(n_falsep/length(RECORDS)) );
fprintf(fileID1, 'False Negative: %d/%d (%f) => Calculated False but Expected True\n', n_falsen, length(RECORDS), 100*(n_falsen/length(RECORDS)) );
fprintf(fileID1, 'Total Files Processed: %d\n', length(RECORDS));
fprintf(fileID1, 'Challenge 2015 Score: %f\n', challenge_score);
fprintf(fileID1,'=========================================================================\n');

fclose(fileID1);
fclose(fileID2);