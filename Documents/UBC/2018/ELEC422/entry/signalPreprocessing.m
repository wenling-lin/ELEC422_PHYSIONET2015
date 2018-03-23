%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ELEC 422 - Physionet 2015 Challenge
% Authors:  Wen-Ling Lin (31120132)
%           Janelle Somerville (52155132)
%           Claire Chang (35340132)
%
% Sources Used: 
% Pre-Processing Filter Specs Used: http://www.cinc.org/archives/2015/pdf/1181.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [filteredSignal] = signalPreprocessing(rawSignal)
rawSignal = rawSignal';
[m,n] = size(rawSignal);
    for c = 1:m
        [C,L]=wavedec(rawSignal(c,:),8,'db4'); % [approx, detail] as per lecture
        [thr,sorh,keepapp]=ddencmp('den','wv',rawSignal(c,:));
        cleanSignal=wdencmp('gbl',C,L,'db4',8,thr,sorh,keepapp);

        % Clean unneeded variables - Prevents the MEMORY Error from occuring
        clear f_y; clear C; clear L;

        % Remove Baseline Wander
        opol = 6; %polynomial degree
       % ecg_signal = ecg1';
        [p,s,mu] = polyfit((1:numel(cleanSignal)),cleanSignal,opol);
        f_y = polyval(p,(1:numel(cleanSignal)),[],mu);
        filteredSignal(:,c) = cleanSignal - f_y;
    end
end

