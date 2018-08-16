function [sync1] = readSync1Data(file)

d = textread(file,'%s');
[date, time, TTL] = textread(file,'%s %s %s');

t_abs = datetime(strcat(date,'_',time),'InputFormat','dd.MM.yyyy_HH:mm:ss.SSSS');
TTL = str2num(cell2mat(TTL));

for l = 2:length(TTL)    
    if(TTL(l-1)==0 & TTL(l)==1)
        sync1 = t_abs(l); % export absolute arrival time of first TTL 1
    end    
end

end