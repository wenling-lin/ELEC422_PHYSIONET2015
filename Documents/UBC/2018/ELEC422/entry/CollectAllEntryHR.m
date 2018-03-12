clc;clear;clear all;

clear all;close all;clc
data_dir=[pwd filesep];

addpath(pwd)

fid=fopen([data_dir 'ALARMS'],'r');
if(fid ~= -1)
    RECLIST=textscan(fid,'%s %s %d','Delimiter',',');
    fclose(fid);
else
    error('Could not open ALARMS.txt. Exiting...')
end

RECORDS=RECLIST{1};
ALARMS=RECLIST{2};
N=length(RECORDS);
results=zeros(N,1);

total_time=0;
for i=1:N
    fname=RECORDS{i};
    tic;
    display(['results(' num2str(i) ')=challenge(''' data_dir fname ''','''  ALARMS{i} ''');'])
    try
        results(i)=challenge([data_dir fname],ALARMS{i});
    catch
        warning(lasterr)
    end
    total_time=total_time+toc;
    if(~mod(i,10))
        fprintf(['---Processed ' num2str(i) ' out of ' num2str(N) ' records.\n'])
    end
end