%% Experiment information

animal = '5627rr';
sessionV = 176; % 176
sessionS = 178; % 178

cutoff = 0.6;


%% Read files and get behavior

bin = 40;

eV = readExperimentData(strcat(animal,'-s',num2str(sessionV),'-exp.txt'));
bV = getBehavior(eV,bin);
bsV = getBehaviorByStimulus(eV);

eS = readExperimentData(strcat(animal,'-s',num2str(sessionS),'-exp.txt'));
bS = getBehavior(eS,bin);
bsS = getBehaviorByStimulus(eS);


%% Find period of proper task performance

startV = nanmax(find(bV.lr_H>=cutoff,1),find(bV.lr_CR>=cutoff,1));
stopV = startV+find(bV.lr_H(startV:end)<cutoff,1)-2; %nanmin(...,find(bV.lr_CR(startV:end)<0.6,1));

startS = nanmax(find(bS.lr_H>=cutoff,1),find(bS.lr_CR>=cutoff,1));
stopS = startS+find(bS.lr_H(startS:end)<cutoff,1)-2; %nanmin(...,find(bS.lr_CR(startS:end)<0.6,1));


%% Calculate metrics for subsets

V = getBehaviorSubset(eV,bV,startV,stopV);
S = getBehaviorSubset(eS,bS,startS,stopS);

V.dp
S.dp
V.n_tot
S.n_tot


%% Plot attention switch

plotAttentionSwitch(bV,bS,bin,startV,stopV,startS,stopS);
%plotAttentionSwitch2(bV,bS,bin,startV,stopV,startS,stopS);

