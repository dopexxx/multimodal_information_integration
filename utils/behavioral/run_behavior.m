% Task: Create compressed version of behavioral data.
% Atm: H:data (518 GB)
% Turned into < 10MB 

root = "H:\data\behavior\";
delim = ";";

files = dir(strcat(root,"*\*\*\*-exp.txt"));

for file_ind = 1:length(files)
    
    file_path = strcat(files(file_ind).folder,"\",files(file_ind).name);
    metastats = readExperimentData(...
        strcat(files(file_ind).folder,'\'), files(file_ind).name ...
    );
    counts = getBehavior(metastats, 40);
    save(strcat(files(file_ind).folder,"\metadata.m"),'metastats');
    save(strcat(files(file_ind).folder,"\response_counts.m"),'counts');

end
