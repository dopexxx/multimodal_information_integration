%% Analyze behavior data

animal = '5627rr'; % '5212r','5627rr'
session = '103';
binSize = 30;

e = readExperimentData(strcat(animal,'-s',session,'-exp.txt')); 
b = getBehavior(e,binSize);

%% Analyze behavior data subset

if isempty(find(b.lr_FA<0.4,1))
    start = 1;
else
    start = find(b.lr_FA<0.4,1);
end

if isempty(find(b.lr_H<0.4,1))
    stop = length(b.lr_H);
else
    stop = find(b.lr_H<0.4,1)+binSize-1;
end

% behavior count

n.beh = e.beh(start:stop);
n.stim = e.stim(start:stop);

n.n_tot = length(n.beh);
n.n_H = 0;
n.n_M = 0;
n.n_FA = 0;
n.n_CR = 0;
for i=1:n.n_tot
    if cell2mat(n.beh(i)) == 'H'
        n.n_H = n.n_H+1;
    elseif cell2mat(n.beh(i)) == 'M'
        n.n_M = n.n_M+1;
    elseif cell2mat(n.beh(i)) == 'FA'
        n.n_FA = n.n_FA+1;    
    elseif cell2mat(n.beh(i)) == 'CR'
        n.n_CR = n.n_CR+1;
    end
end

% global rates

if n.n_H+n.n_M==0
    n.n_go = 1;
else
    n.n_go = n.n_H+n.n_M;
end
if n.n_FA+n.n_CR==0
    n.n_nogo = 1;
else
    n.n_nogo = n.n_FA+n.n_CR;
end
if n.n_H+n.n_FA==0
    n.n_lick = 1;
else
    n.n_lick = n.n_H+n.n_FA;
end
if n.n_M+n.n_CR==0
    n.n_nolick = 1;
else
    n.n_nolick = n.n_M+n.n_CR;
end

n.r_H = n.n_H/n.n_go;  
n.r_M = n.n_M/n.n_go;  
n.r_FA = n.n_FA/n.n_nogo;  
n.r_CR = n.n_CR/n.n_nogo;

if n.r_H == 1
    n.r_H = 0.99;
    n.r_M = 0.11;
    disp('Values changed for calculation of d prime')
end
[n.dp,n.c] = dprime_simple(n.r_H,n.r_FA);


%% Plot behavior data

figure

hold on
x = [30:length(b.lr_H)+29];
plot(x,b.lr_H,'g')
plot(x,b.lr_CR,'c')
line([start+29 start+29], get(gca,'ylim'),'Color','r');
line([stop-binSize+1+29 stop-binSize+1+29], get(gca,'ylim'),'Color','r');
xlim([30 length(b.lr_H)+29])
ylim([0 1])
xlabel('trial')
ylabel('rate')
legend('H','CR')
hold off

