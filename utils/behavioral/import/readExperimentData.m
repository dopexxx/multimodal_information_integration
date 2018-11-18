function [e] = readExperimentData(path,filename)

% deprecated. Use strcat, which works with both, char and string.s - Jannis
%full=[path,filename] 
full = strcat(path,filename)
    
fid = fopen(full); % e.g. '5627rr-s3-exp.txt'
C = textscan(fid, '%s','delimiter', '\t');
fclose(fid);

[truefalse, idx_1] = ismember('EXPERIMENT INFORMATION', C{1,1});
[truefalse, idx_2] = ismember('REWARD AND PUNISHMENT', C{1,1});
[truefalse, idx_3] = ismember('EXPERIMENT STRUCTURE', C{1,1});
[truefalse, idx_4] = ismember('TRIAL TYPES', C{1,1});
[truefalse, idx_5] = ismember('TRIAL SEQUENCE', C{1,1});
[truefalse, idx_6] = ismember('BEHAVIOR SEQUENCE', C{1,1});
[truefalse, idx_7] = ismember('TRIAL SELECTION SEQUENCE', C{1,1});

%e.t_abs = datetime(strcat(C{1,1}(idx_1+1),C{1,1}(idx_1+2)),'InputFormat','dd.MM.yyyy HH:mm:ss.SSSS');
e.animal = cell2mat(C{1,1}(idx_1+4));
e.session = str2num(cell2mat(C{1,1}(idx_1+6)));
e.n_trials = str2num(cell2mat(C{1,1}(idx_1+8)));
e.freeze = str2num(cell2mat(C{1,1}(idx_1+10)));

e.t_open = str2num(cell2mat(C{1,1}(idx_2+2)));
e.cue = str2num(cell2mat(C{1,1}(idx_2+4)));
e.autoReward = str2num(cell2mat(C{1,1}(idx_2+6)));
e.punishFalseAlarm = str2num(cell2mat(C{1,1}(idx_2+8)));
e.punishEarlyLicks = str2num(cell2mat(C{1,1}(idx_2+10)));

e.t_pre = str2num(cell2mat(C{1,1}(idx_3+2)));
e.t_in = str2num(cell2mat(C{1,1}(idx_3+4)));
e.t_stim = str2num(cell2mat(C{1,1}(idx_3+6)));
e.t_out = str2num(cell2mat(C{1,1}(idx_3+8)));
e.t_post = str2num(cell2mat(C{1,1}(idx_3+10)));
e.t_rew = str2num(cell2mat(C{1,1}(idx_3+12)));
e.t_inter = str2num(cell2mat(C{1,1}(idx_3+14)));

e.f_equal = str2num(cell2mat(C{1,1}(idx_4+2)));
e.N = str2num(cell2mat(C{1,1}(idx_4+4)));
e.N_ratio = str2num(cell2mat(C{1,1}(idx_4+5)));
e.V = str2num(cell2mat(C{1,1}(idx_4+7)));
e.V_ratio = str2num(cell2mat(C{1,1}(idx_4+8)));
e.S = str2num(cell2mat(C{1,1}(idx_4+10)));
e.S_ratio = str2num(cell2mat(C{1,1}(idx_4+11)));
e.VS = str2num(cell2mat(C{1,1}(idx_4+13)));
e.VS_ratio = str2num(cell2mat(C{1,1}(idx_4+14)));
e.Sd = str2num(cell2mat(C{1,1}(idx_4+16)));
e.Sd_ratio = str2num(cell2mat(C{1,1}(idx_4+17)));
e.learnCorrectRejection = str2num(cell2mat(C{1,1}(idx_4+19)));
e.enforceCorrect = str2num(cell2mat(C{1,1}(idx_4+21)));

e.stim = C{1,1}(idx_5+1:idx_6-1);
e.beh = C{1,1}(idx_6+1:idx_7-1);
e.sel = C{1,1}(idx_7+1:end);

end