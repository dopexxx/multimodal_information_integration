%% Calculations

clear

animal = '5212r'; % '5212r','5627rr'
sessions = 230;

% [~,list]=system('dir /S *.txt');

a.lr_H = [];
a.lr_CR = [];
a.l_dp = [];
a.l_c = [];

for session = 1:sessions
    
    % preparations
    sessionList{session} = strcat('s',num2str(session));
    files.e.(sessionList{session}) = strcat(animal,'-s',num2str(session),'-exp.txt');
    
    % read files
    e.(sessionList{session}) = readExperimentData(files.e.(sessionList{session}));
    
    % get behavior
    b.(sessionList{session}) = getBehavior(e.(sessionList{session}));
    
    % save H rate, CR rate and d'
    a.r_H(session) = b.(sessionList{session}).r_H;
    a.r_CR(session) = b.(sessionList{session}).r_CR;
    a.dp(session) = b.(sessionList{session}).dp;
    a.c(session) = b.(sessionList{session}).c;
    a.lr_H = vertcat(a.lr_H, b.(sessionList{session}).lr_H);
    a.lr_CR = vertcat(a.lr_CR, b.(sessionList{session}).lr_CR);
    a.l_dp = vertcat(a.l_dp, b.(sessionList{session}).l_dp);
    a.l_c = vertcat(a.l_c, b.(sessionList{session}).l_c);
end


%% Plots

figure

subplot(2,2,1)
hold on
plot(a.r_H,'g')
plot(a.r_CR,'c')
ylim([0 1])
title(animal)
xlabel('session')
ylabel('rate')
legend('H','CR')
hold off

subplot(2,2,3)
hold on
plot(a.dp)
plot(a.c)
title(animal)
xlabel('session')
ylabel('standard deviations')
legend('sensitivity: d prime','response bias: c')
hold off

subplot(2,2,2)
hold on
plot(a.lr_H,'g')
plot(a.lr_CR,'c')
ylim([0 1])
title(animal)
xlabel('trial')
ylabel('local rate')
legend('H','CR')
border = 0;
for session = 1:sessions
    border = border + length(b.(sessionList{session}).lr_H);
    line([border border], get(gca,'ylim'),'Color','k');
end
hold off

subplot(2,2,4)
hold on
plot(a.l_dp)
plot(a.l_c)
title(animal)
xlabel('trial')
ylabel('standard deviations')
legend('sensitivity: local d prime', 'response bias: local c')
border = 0;
for session = 1:sessions
    border = border + length(b.(sessionList{session}).l_dp);
    line([border border], get(gca,'ylim'),'Color','k');
end
hold off