%% Read data

clear

e = readExperimentData('Fridolin-s3-exp.txt');
t = readTrialData('Fridolin-s3-t1.txt');
l = readLickData('Fridolin-s3-t1-l.txt');


%% Align lick and trial data

if t.t_abs(1) < l.t_abs(1)
    a.t0 = t.t_abs(1);
else
    disp('ERROR: Lick data start earlier than trial data');
end
a.t_diff = duration(l.t_abs(1)-a.t0,'Format','s');

a.t_lick = l.t_rel + a.t_diff;


%% Find beginning of new state

[~, a.idx_pre] = ismember('pre', t.state);
[~, a.idx_in] = ismember('in', t.state);
[~, a.idx_stim] = ismember('stim', t.state);
%[~, a.idx_out] = ismember('out', t.state);
[~, a.idx_post] = ismember('post', t.state);
[~, a.idx_rew] = ismember('rew', t.state);
[~, a.idx_inter] = ismember('inter', t.state(2:end));

a.t_pre = t.t_rel(a.idx_pre);
a.t_in = t.t_rel(a.idx_in);
a.t_stim = t.t_rel(a.idx_stim);
%a.t_out = t.t_rel(a.idx_out);
a.t_post = t.t_rel(a.idx_post);
a.t_rew = t.t_rel(a.idx_rew);
a.t_inter = t.t_rel(a.idx_inter);


%% Calculate real state durations

a.d_pre = a.t_in - a.t_pre;
a.d_in = a.t_stim - a.t_in;
%a.d_stim = a.t_out - a.t_stim;
a.d_stim = a.t_post - a.t_stim; % instead
%a.d_out = a.t_post - a.t_out;
a.d_post = a.t_rew - a.t_post;
a.d_rew = a.t_inter - a.t_rew;

% calculate duration since last time point -> median ~ 0.064 s
a.ltp = double.empty(length(t.t_rel)-1,0);
for i = 2:length(t.t_rel)
    a.ltp(i-1) = seconds(t.t_rel(i) - t.t_rel(i-1));
end


%% Plot

figure

hold on
plot(a.t_lick,l.lick)
title(strcat(e.animal,'-s',num2str(e.session)))
xlabel('time [s]')
ylabel('lick behavior')
line([seconds(a.t_stim),seconds(a.t_stim)], get(gca,'ylim'),'Color','r','LineWidth',3);
%line([seconds(a.t_out),seconds(a.t_out)], get(gca,'ylim'),'Color','r','LineWidth',3);
line([seconds(a.t_rew),seconds(a.t_rew)], get(gca,'ylim'),'Color','b','LineWidth',3);
line([seconds(a.t_inter),seconds(a.t_inter)], get(gca,'ylim'),'Color','b','LineWidth',3);
hold off

