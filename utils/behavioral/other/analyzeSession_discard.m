%% Read file

animal = '5212r'; % '5212r','5627rr'
session = '201';

e = readExperimentData(strcat(animal,'-s',session,'-exp.txt'));


%% Get behavior
   
b = getBehavior(e,40);

disp(strcat('H:',num2str(b.r_H)));
disp(strcat('CR:',num2str(b.r_CR)));
disp(strcat('d prime:',num2str(b.dp)));
disp(strcat('c:',num2str(b.c)));


%% Cut beginning and end

if isempty(find(b.lr_FA<0.4,1))
    start = 1;
else
    start = find(b.lr_FA<0.4,1)
end
start = find(b.lr_FA<0.4,1)

if isempty(find(b.lr_H<0.4,1))
    stop = length(b.lr_H);
else
    stop = find(b.lr_H<0.4,1)
end


%% Plots

figure

subplot(2,1,1)
hold on
plot(b.lr_H,'g')
plot(b.lr_CR,'c')

xlim([0 length(b.lr_H)])
ylim([0 1])
title(animal)
xlabel('trial')
ylabel('rate')
legend('H','CR')
line([start start], get(gca,'ylim'),'Color','k');
line([stop stop], get(gca,'ylim'),'Color','k');
hold off

subplot(2,1,2)
hold on
plot(b.l_dp)
plot(b.l_c)
xlim([0 length(b.lr_H)])
title(animal)
xlabel('trial')
ylabel('standard deviations')
legend('sensitivity: d prime', 'response bias: c')
hold off