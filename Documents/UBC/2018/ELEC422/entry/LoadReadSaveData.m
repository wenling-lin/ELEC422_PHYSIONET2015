clear; clc; clear all;

% ELEC 422 - Biosignals
% Purpose: Save this script in the same folder as where the data is saved and
% it will filter, load and save plots into a specified directory.

% Prompt user to filter data set by the first letter of the files
prompt = 'Enter first letter of the type of files to process (Valid Inputs: v, t, f, b, a)\n';
usr_input = input(prompt, 's');
fprintf('\n');

% Prompt user to select output directory
prompt = 'Enter output directory to save plots\n';
output_dir = input(prompt, 's');
fprintf('\n');
current_dir = pwd; % Save current directory

%Create Directory if directory does not exist
if( exist(output_dir)== 0 )
    fprintf('Directory does not exist, creating directory...\n\n');
    mkdir (output_dir);
end

% Search for files with inputted prefix
files_with_prefix = strcat(usr_input, '*', '.mat');
filtered_files=dir(files_with_prefix);

%Load the data and save files
num_files=size(filtered_files);
num_files=num_files(1);

for i = 1:num_files
    cd (current_dir)
    % Load data
    data=load(filtered_files(i).name);
    msg = strcat({'Loading File: '},{filtered_files(i).name});
    disp(msg);
    data=data.val;
    
    % Check how many biosignals in loaded file
    num_biosignals=size(data(:,1));
    num_biosignals=num_biosignals(1);
    
    %figure('units','normalized','outerposition',[0 0 1 1])
    h=figure('Name', filtered_files(i).name,'units','normalized','outerposition',[0 0 1 1]);
    % PPG Signal
    subplot(4,1,1)
    plot(data(1,:))
    title('PPG')
    
    %ABP
    subplot(4,1,2)
    plot(data(2,:))
    title('ABP')

    %ECG
    subplot(4,1,3)
    plot(data(3,:))
    title('ECG')

    %ECGAVR -- (Plots fourth plot if fourth signal exists)
    if( num_biosignals > 3 )
        subplot(4,1,4)
        plot(data(4,:))
        title('ECGAVR')
    end
    
    %Save figure if figures do not exist
    cd (output_dir);
    filename_fig = strcat(filtered_files(i).name,'.fig');
    filename_png = strcat(filtered_files(i).name,'.png');
    if( exist(filename_fig) == 0 )
        msg = strcat({'Saving Plot: '},{filtered_files(i).name});
        disp(msg);
        savefig(filename_fig);
        saveas(h,filename_png);
    end
    %Close plot
    close(h)
    
end

cd (current_dir)





