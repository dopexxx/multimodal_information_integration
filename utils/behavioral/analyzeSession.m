%% Select data

animal = '5212r'; % e.g.'5627rr'
session = 227; % e.g. 347

bin = 40; %smoothing window for local performance

%% Behavioral analysis

e = readExperimentData(strcat(animal,'-s',num2str(session),'-exp.txt'));

b = getBehavior(e,bin);

disp(strcat('n:',num2str(b.n_tot)));
disp(strcat('H:',num2str(b.r_H)));
disp(strcat('CR:',num2str(b.r_CR)));
disp(strcat('d prime:',num2str(b.dp)));
disp(strcat('EL:',num2str(b.r_EL)));
disp(strcat('correct:',num2str(b.r_correct)));

%% Plot

figure
hold on
plot(b.lr_H,'g')
plot(b.lr_CR,'c')
plot(b.lr_EL,'y')
xlim([0 length(b.lr_H)])
ylim([0 1])
title(num2str(session))
xlabel('trial')
ylabel('rate')
legend('H','CR','EL')
hold off