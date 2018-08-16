%% Parameters

animal = '5212r';
sessionV = 339;

bin = 40;
cutoff = 0.6;

%% Processing

eV = readExperimentData(strcat(animal,'-s',num2str(sessionV),'-exp.txt'));
bV = getBehavior(eV,bin);
bsV = getBehaviorByStimulus(eV);

startV = nanmax(find(bV.lr_H>=cutoff,1),find(bV.lr_CR>=cutoff,1));
stopV = startV+find(bV.lr_H(startV:end)<cutoff,1)-2; %nanmin(...,find(bV.lr_CR(startV:end)<0.6,1));
V = getBehaviorSubset(eV,bV,startV,stopV);

V.dp
V.n_tot

%% Plot

figure
hold on
xV = [bin:length(bV.lr_H)+(bin-1)]';
plot(xV,bV.lr_H,'g')
plot(xV,bV.lr_CR,'c')
plot(xV,bV.lr_EL,'y')
xlim([bin,length(bV.lr_H)+(bin-1)])
ylim([0,1])
line(get(gca,'xlim'),[0.6,0.6],'Color','k','LineStyle','--');
line([40,40],get(gca,'ylim'),'Color','r');
line([271,271],get(gca,'ylim'),'Color','r');
%line([167,167],get(gca,'ylim'),'Color','r');
%line([217,217],get(gca,'ylim'),'Color','r');
%line([startV+(bin-1),startV+(bin-1)],get(gca,'ylim'),'Color','r');
%line([stopV+(bin-1),stopV+(bin-1)],get(gca,'ylim'),'Color','r');
title('202')
%a = area([0,86],[1,1]);
%a.FaceAlpha = 0.2;
xlabel('trial')
ylabel('rate')
legend('H','CR','EL')
hold off