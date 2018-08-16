function sync1 = loadSync1(date, animal, session, trials)

% get path
if ispc
    % path = strcat('D:/data/behavior/',date,'/',animal,'/',num2str(session),'/');
    disp('SPECIFY PATH OF HD IN loadSync1')
elseif ismac
    path = strcat('/Volumes/Moritz_body/data/behavior/',date,'/',animal,'/',num2str(session),'/');
end

for trial = 1:trials
    % load data from file
    [date1, time, TTL] = textread(strcat(path,animal,'-s',num2str(session),'-t',num2str(trial),'-s.txt'),'%s %s %s');

    t_abs = datetime(strcat(date1,'_',time),'InputFormat','dd.MM.yyyy_HH:mm:ss.SSSS');
    TTL = str2num(cell2mat(TTL));

    % find absolute time of TTL onsets for each trial
    for l = 2:length(TTL)    
        if(TTL(l-1)==0 && TTL(l)==1)
            sync1_temp{trial} = t_abs(l);
        end    
    end
    disp(['Successfully loaded sync of trial ', num2str(trial), '.'])
end

sync1 = [sync1_temp{:}]';

end


