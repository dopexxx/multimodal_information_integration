function [rel,stim,beh,sync1] = loadExperiment(date, animal, session)

% get path
if ispc
    path = strcat('G:/data/behavior/',date,'/',animal,'/',num2str(session),'/');
elseif ismac
    %path = strcat('/Users/jannisborn/Desktop/HIFO/behavior/',date,'/',animal,'/',num2str(session),'/');
    path = '/Users/jannisborn/Desktop/HIFO/227/';

end

fid = fopen(strcat(path,animal,'-s',num2str(session),'-exp.txt'));
C = textscan(fid, '%s','delimiter', '\t');
fclose(fid);

[~,idx_1] = ismember('EXPERIMENT INFORMATION', C{1,1});
[~,idx_5] = ismember('TRIAL SEQUENCE', C{1,1});
[~,idx_6] = ismember('BEHAVIOR SEQUENCE', C{1,1});
[~,idx_7] = ismember('TRIAL SELECTION SEQUENCE', C{1,1});
[~,idx_sync1] = ismember('SYNC  SEQUENCE', C{1,1});
%rel = C{1,1}(...);
data_stim = C{1,1}(idx_5+1:idx_6-1);
data_beh = C{1,1}(idx_6+1:idx_7-1);
data_sync1 = C{1,1}(idx_sync1+1:end);

for trial = 1:str2num(cell2mat(C{1,1}(idx_1+8)));
    trialString = strcat('t',num2str(trial));
    
    stim.(trialString) = data_stim(trial);
    beh.(trialString) = data_beh(trial);
    sync1_temp.(trialString) = datetime(data_sync1(trial),'InputFormat','dd.MM.yyyy HH:mm:ss.SSSS');
    sync1(trial) = sync1_temp.(trialString);
end
sync1 = sync1';
%sync1 = 'error'
    
% modality relevance

[~,idx_rel] = ismember('MODALITY RELEVANCE', C{1,1});
if str2num(cell2mat(C{1,1}(idx_rel+2)))==1 && str2num(cell2mat(C{1,1}(idx_rel+4)))==0
    rel = 'V';
elseif str2num(cell2mat(C{1,1}(idx_rel+2)))==0 && str2num(cell2mat(C{1,1}(idx_rel+4)))==1
    rel = 'S';
elseif str2num(cell2mat(C{1,1}(idx_rel+2)))==0 && str2num(cell2mat(C{1,1}(idx_rel+4)))==0
    rel = 'V'; % of course not true
end

end