%% Plot results of motion analysis
load('motion_analysis.mat')
close all
figure
for k = 1:4
    plot(1:length(motion.Results),squeeze(motion.Results(1,k,:))/max(motion.Results(1,k,:)));
    hold on;
end
legend(motion.algos)
motion.algos
%% Show correlation coefficient of methods
for k = 1:length(motion.ROIs)
    strcat("CorrCoef matrix of ROI ", motion.ROIs{k}) 
    corrcoef(squeeze(motion.Results(k,:,:))')
end




%% Function testing
close all;clc;
path = pwd;
test_motion_algo(pwd,13,1)


%%
figure
plot(squeeze(motion.Results(1,2,:)))
