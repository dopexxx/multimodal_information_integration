%% Read files

clear

animal = '5627rr'; % '5212r','5627rr'
sessions = [8, 9];

sessionList = num2str([]);

for session = 1:length(sessions)
    
    % preparations
    sessionList = [sessionList, strcat('s',num2str(sessions(session)))];
    files.e.(sessionList{session}) = strcat(animal,'-s',num2str(session),'-exp.txt');
    
    % read files
    e.(sessionList{session}) = readExperimentData(files.e.(sessionList{session}));   
end

