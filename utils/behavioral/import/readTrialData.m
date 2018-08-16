function [t] = readTrialData(filename)

fid = fopen(filename); % e.g. '5212r-s4-t1.txt'
C = textscan(fid, '%s %s %s %s %s %s %s %s %s','delimiter', '\t');
fclose(fid);

t.t_abs = datetime(strcat(C{1,1},C{1,2}),'InputFormat','dd.MM.yyyy HH:mm:ss.SSSS');
t.t_rel = duration(t.t_abs - t.t_abs(1), 'Format', 's');
t.trial = str2num(cell2mat(C{1,3}(1)));
t.type = cell2mat(C{1,4}(1));
t.enum = C{1,5};
t.state = C{1,6};
t.V = str2num(cell2mat(C{1,7}));
t.S = str2num(cell2mat(C{1,8}));
t.stage = str2double(C{1,9});

end