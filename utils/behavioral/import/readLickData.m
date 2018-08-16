function [l] = readLickData(filename)

[date, time, lick] = textread(filename,'%s %s %s'); % e.g. 'Fridolin-s21-t1-l.txt'

l.t_abs = datetime(strcat(date,'_',time),'InputFormat','dd.MM.yyyy_HH:mm:ss.SSSS');
l.t_rel = duration(l.t_abs - l.t_abs(1), 'Format', 's');
l.lick = str2num(cell2mat(lick));

end