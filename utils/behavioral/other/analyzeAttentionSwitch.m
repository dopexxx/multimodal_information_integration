%% Experiment information

animal = '5627rr';
sessionV = 71;
sessionS = 73;
% 12.07.: 65&67
% 13.07.: 71&73
% 14.07.: 75&77
% 14.07.: 79&81
% 15.07.: 83&85
% 17.07.: 92&94 okay
% best: 176,178

%% Read files and get behavior

bin = 40;

eV = readExperimentData(strcat(animal,'-s',num2str(sessionV),'-exp.txt'));
bV = getBehavior(eV,bin);
bsV = getBehaviorByStimulus(eV);

eS = readExperimentData(strcat(animal,'-s',num2str(sessionS),'-exp.txt'));
bS = getBehavior(eS,bin);
bsS = getBehaviorByStimulus(eS);

%% Compare relevance conditions - numbers

fprintf(['\nTotal number of trials: ',num2str(eV.n_trials),' Vrel vs. ',num2str(eS.n_trials),' Srel.\n\n'])

fprintf(['Hit rate: ',num2str(bV.r_H),' Vrel vs. ',num2str(bS.r_H),' Srel.\n'])
fprintf(['Correct rejection rate: ',num2str(bV.r_CR),' Vrel vs. ',num2str(bS.r_CR),' Srel.\n'])
fprintf(['Early lick rate: ',num2str(bV.r_EL),' Vrel vs. ',num2str(bS.r_EL),' Srel.\n'])
fprintf(['Correct trials rate: ',num2str(bV.r_correct),' Vrel vs. ',num2str(bS.r_correct),' Srel.\n'])
fprintf(['d prime: ',num2str(bV.dp),' Vrel vs. ',num2str(bS.dp),' Srel.\n\n'])

fprintf(['Correct trials rate for relevant stimulus alone: ',num2str(bsV.r_Vc),' Vrel vs. ',num2str(bsS.r_Sc),' Srel.\n'])
fprintf(['Early lick rate for relevant stimulus alone: ',num2str(bsV.r_VEL),' Vrel vs. ',num2str(bsS.r_SEL),' Srel.\n\n'])

fprintf(['Correct trials rate for irrelevant stimulus alone: ',num2str(bsV.r_Sc),' Vrel vs. ',num2str(bsS.r_Vc),' Srel.\n'])
fprintf(['Early lick rate for irrelevant stimulus alone: ',num2str(bsV.r_SEL),' Vrel vs. ',num2str(bsS.r_VEL),' Srel.\n\n'])

%% Compare relevance conditions - plots
figure

subplot(2,1,1)
hold on
plot(bV.lr_H,'g')
plot(bV.lr_CR,'c')
plot(bV.lr_EL,'r')
xlim([0 length(bV.lr_H)])
ylim([0 1])
title('V relevant')
xlabel('trial')
ylabel('rate')
legend('H','CR','EL')
hold off

subplot(2,1,2)
hold on
plot(bS.lr_H,'g')
plot(bS.lr_CR,'c')
plot(bS.lr_EL,'r')
xlim([0 length(bS.lr_H)])
ylim([0 1])
title('S relevant')
xlabel('trial')
ylabel('rate')
legend('H','CR','EL')
hold off
