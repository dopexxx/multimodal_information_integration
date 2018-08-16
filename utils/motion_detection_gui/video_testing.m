
%%
filename = '/Volumes/Moritz_body/data/body/17-10-03/5627rr/208/good.avi';

vid = VideoReader(filename);
chunked = VideoWriter('test','mpeg-4');
open(chunked);


for r = 1:50
    new = read(vid,r);
    imshow(new)
    writeVideo(chunked,new);
end
close(chunked);
close all;

%%
old = rgb2gray(read(vid,3)); % first frame

for ind = 4:4
    new = rgb2gray(read(vid,ind));
end
tic
opticFlow = opticalFlowLK('NoiseThreshold',0.009); % optical flow estimator specifications 
toc
flow = estimateFlow(opticFlow,new);
motion = flow.Magnitude;
toc
%%

a = ['Here you can select which motion detection algorithms should be applied.', ...
    newline,'Per default, all will be applied. Untick to turn them off and speed up processing. Confirm with okay.'];
            

r = zeros(size(motion));
%%
% Read and display a MP4 file
filename = '/Users/jannisborn/Desktop/HIFO/body/17-06-09/5627rr/101/5627rr-s101-body_c.mp4';

v = VideoWriter('test/test.mp4','mpeg-4');
vv = VideoReader(filename);
open(v)
% Create video
old = read(vv,1);

move = zeros([stopp-start+1,size(old)]);

start = 2;
stopp = 100;
count = 0;

for ind = start:stopp
    count = count+1;
    new = read(vv,ind);
    diff = imabsdiff(old,new);
    %diff = new;
    move(count,:,:,:) = diff;
    writeVideo(v,diff);
    old = new;
    imagesc(diff);
    pause(0.2)
end
close(v)

42

%%
vv = VideoReader(filename);
opticFlow = opticalFlowLK('NoiseThreshold',0.009);
%opticFlow = opticalFlowHS;
start = 2;
stopp = 200;
count = 0;

for ind = start:stopp
    %frameRGB = readFrame(vv);
    frameRGB = read(vv,ind);
    frameGray = rgb2gray(frameRGB);
    flow = estimateFlow(opticFlow,frameGray); 

%     imshow(frameRGB) 
%     hold on
%     plot(flow,'DecimationFactor',[5 5],'ScaleFactor',10)
%     hold off
    imagesc(flow.Magnitude)
    pause(0.2)
end

%%
nFrames = ceil(video.FrameRate*video.Duration)
vidHeight = video.Height
vidWidth = video.Width
fps = round(video.FrameRate)
%%
initial = 1;
final = 1000
frameNum = final-initial;
mov(1:frameNum) = struct ('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'), 'colormap', []);
for i = 1:frameNum
    mov(i).cdata= read(video,((initial-1)+i));
    image(mov(i).cdata)
end
%%
numFrames = ceil(v.FrameRate*v.Duration)
%%
currAxes = axes;
k=0;
while hasFrame(v)
    k=k+1;
    if mod(k,100)==0
        k
    end
    tic 
    vidFrame = readFrame(v,k);
    toc
    image(vidFrame, 'Parent', currAxes);
    currAxes.Visible = 'off';
    pause(1/v.FrameRate);
end


%%
tic
filename = '/Users/jannisborn/Desktop/HIFO/Videos/5212r-s101-body_c.mp4';
video = VideoReader(filename);

% Compute metadata
nFrames = ceil(video.FrameRate*video.Duration);
vidHeight = video.Height;
vidWidth = video.Width;
fps = round(video.FrameRate);

% Temporal cropping (which frames to present)
initial = 1;
final = 2 * fps;
mov(1:final) = struct ('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'), 'colormap', []);
toc

%%
for i = 1:final
    mov(i).cdata= read(video,((initial-1)+i));
    imshow(mov(i).cdata)
end


%%
close all
vidSrc = vision.VideoFileReader(filename,'ImageColorSpace','RGB');
img = step(vidSrc);
imshow(img);
% drawnow();
while ~isempty(img)
    imshow(step(vidSrc));
%     drawnow();
end

%%
close all; clc;
video = vision.VideoFileReader(filename,'ImageColorSpace','RGB');
 
first_frame = step(video);
imshow(first_frame);

final = 10 * 20; % 20 frames per second in all videos

for ind = 1:final
    imshow(step(video))
    drawnow
end

%%
path = '/Users/jannisborn/Desktop/HIFO/body/17-06-09/5212r/101/5212r-s101-body_c.mp4';

% Detection requests that 
%   1) the date follows the '/' after the only 'body' string in all but the
%           last 10 chars of the string
%   2) the date has 7! chars
%   3) the file name starts with the animal name followed by '-'
%   4) the session name follows the animal name and is prefixed and
%       suffixed by '-' 

date_onset = strfind(path(1:length(path)-10),'body')+5;
date = path(date_onset:date_onset+7); % date has fixed length

animal_on = strfind(path,'/'); animal_on=animal_on(end)+1;
% animal has not fixed length, animal_off is relative to animal_on
animal_off = strfind(path(animal_on:end),'-'); 
animal_off = animal_off(1) + animal_on - 2;
animal = path(animal_on:animal_off);

% session follows straight after the animal, separated by -
session_on = animal_off + 2;
session_off = strfind(path(session_on:end),'-'); 
session_off = session_on + session_off(1) - 2; % in case > 1 '-' in filename
session = path(session_on:session_off);


%%

cd '/Users/jannisborn/Desktop/HIFO/body';
[baseName, folder] = uigetfile('.mp4')
fullFileName = fullfile(folder, baseName)

%%
filename = '5212r-s227-t1.txt';
fid = fopen(filename); % e.g. '5627rr-s3-exp.txt'
C = textscan(fid, '%s','delimiter', '\t');
fclose(fid);

idx = find(contains(C{1},'stim')) - 4;
day_time = [C{1}{idx-1},C{1}{idx}];
a = datetime(day_time,'Format','dd.MM.yyyy HH:mm:ss.SSSS')


%% duration
clc
trial = 4;
animal = '5212r'; % e.g.'5627rr'
session = '227'; % e.g. 347
fid = fopen(strcat(animal,'-s',session,'-t',num2str(trial),'.txt'));
C = textscan(fid, '%s','delimiter', '\t');

on = datetime([C{1}{1},C{1}{2}],'Format','dd.MM.yyyy HH:mm:ss.SSSS')
off = datetime([C{1}{end-8},C{1}{end-7}],'Format','dd.MM.yyyy HH:mm:ss.SSSS')



%%
clc; close all;
Z = peaks;
surf(Z); 
axis tight manual 

v = VideoWriter('peaks','mpeg-4');
open(v);
for k = 1:20 
   surf(sin(2*pi*k/20)*Z,Z)
   frame = getframe(gcf);
   frame
   writeVideo(v,frame);
end

close(v);
set(gca,'nextplot','replacechildren'); 
clc

%%
f = uifigure;
d = uiprogressdlg(f,'Title','Please Wait','Message','Prepare chunking');

%%
f = uifigure;
    d = uiprogressdlg(f,'Title','Please Wait',...
        'Message','Opening the application');
        pause(.5)
%%
clc
filepath = '/Users/jannisborn/Documents/MATLAB/myBinomTest.m';
target = '/Users/jannisborn/Documents/MATLAB/ledalab-master/';

movefile(filepath, target)

%%
%file = '/Users/jannisborn/Desktop/HIFO/body/17-06-09/5212r/101/5212r-s101-body_c.mp4';
%a=VideoReader(file);
clc
test1 = VideoReader('rhinos.avi')

file = '/Users/jannisborn/Desktop/HIFO/body/17-09-22/5212r/201/good_c.mp4';
test2=VideoReader(file)


%%
close all;clc;
fig = uifigure;
message = ['Please select which of the options from the dropdown you used as a reference. Preferred option', ...
    ' is visual reflection onset. If you choose the stage, make sure the stage has fully moved, i.e. you',...
    ' saw either a S or a V+S trial but not a N or a V trial. This method FAILS otherwise.',...
    ' Confirming will chop the entire video into chunks for every trial.'];
dd = uidropdown(fig, 'Items',{'Onset of visual reflection', 'Offset of stage drive-in'},...
    'ItemsData', {'refl','stage'});
dd.Position = [100 100 170 22];

selection = uiconfirm(fig,message,'Confirm choice','Options',{'OK', 'Cancel'}, 'Defaultoption', 2, ...
'Canceloption',2,'Icon','Warning');


%%
path = '/Users/jannisborn/Desktop/HIFO/behavior/17-06-12/5212r/103/5212r-s103-t3.txt';
fid = fopen(path); % e.g. '5627rr-s3-t2.txt'
C = textscan(fid, '%s','delimiter', '\t');
fclose(fid);
%%

fig = uifigure;
message = sprintf(['You will now draw lines around the regions of interest (ROI) for motion analysis. There ', ...
    'will be a ROI for 1) the snout, 2) the left front paw, 3) the right front paw and 4) the back  ']);

uiwait(msgbox(message));

%%
close all;
fig = uifigure;
handle = gca;

figure
imshow('pout.tif')


%%
close all; clc;
fig = figure('units','normalized','outerposition',[0 0 1 1]);
ax = axes(fig);
imshow('cameraman.tif','parent',ax,'InitialMagnification','fit');
imfreehand(ax)

%%
close all;
img = imread('cameraman.tif');
for k =1:length(a.Yin)
    img(round(a.Yin(k)), round(a.Xin(k))) = 255;
end
imshow(img)


%%
 a = struct();
a.Names = {'Head', 'Left front paw', 'Right front paw', 'Back'};

for k = 1:length(a.Names)
    a.Values{k} = 0;
end
% Just allocation, will be replaced
%a.Values = {0,0,0,0};
%a.raw = {0,0,0,0};
a
            
%%
fig = uifigure;
cb_ad = uicheckbox(fig,'Text','Absolute difference','Value',1,'Position',[150 175 150 15]);
cb_lk = uicheckbox(fig,'Text','Lucas Kanade','Value',0,'Position',[150 150 102 15]);
cb_hs = uicheckbox(fig,'Text','Horn Schunck','Value',0,'Position',[150 125 102 15]);
label1 = uibutton(fig,...
    'Position',[100 50 50 25],...
    'Text','DONE!');

cb_hs



























%% IMABSDIFF
% Read and display a MP4 file
filename = '/Users/jannisborn/Desktop/HIFO/body/17-06-09/5627rr/101/5627rr-s101-body_c.mp4';

v = VideoWriter('test/absdiff.mp4','mpeg-4');
vv = VideoReader(filename);
open(v)
% Create video
old = read(vv,1);

start = 2;
stopp = 100;
count = 0;

for ind = start:stopp
    new = read(vv,ind);
    diff = imabsdiff(old,new);
    writeVideo(v,diff);
    old = new;
    imagesc(diff);
    pause(0.05)
end
close(v)
close all
%% KANADE
filename = '/Users/jannisborn/Desktop/HIFO/body/17-06-09/5627rr/101/5627rr-s101-body_c.mp4';
v = VideoWriter('test/lucasKanade.mp4','mpeg-4');
vv = VideoReader(filename);
open(v)
% Create video
old = rgb2gray(read(vv,1));

start = 2;
stopp = 100;
count = 0;
opticFlow = opticalFlowLK('NoiseThreshold',0.009); % optical flow estimator specifications 

for ind = start:stopp
    new = rgb2gray(read(vv,ind));
    flow = estimateFlow(opticFlow,new);
    motion = flow.Magnitude;
    writeVideo(v,motion/max(motion(:)));
    old = new;
    imagesc(motion);
    pause(0.05)
end
close(v)
close all

%% hornschunk
v = VideoWriter('test/hornschunk.mp4','mpeg-4');
vv = VideoReader(filename);
open(v)
% Create video
old = rgb2gray(read(vv,1));

start = 2;
stopp = 100;
count = 0;
opticFlow = opticalFlowHS;

for ind = start:stopp
    new = rgb2gray(read(vv,ind));
    flow = estimateFlow(opticFlow,new);
    motion = flow.Magnitude;
    writeVideo(v,motion/max(motion(:)));
    old = new;
    imagesc(motion);
    pause(0.05)
end
close(v)

%%
h = animatedline;
axis([0 4*pi -1 1])
x = linspace(0,4*pi,2000);

for k = 1:length(x)
    tic
    y = sin(x(k));
    addpoints(h,x(k),y);
    drawnow
    toc
end
%%
path = '/Volumes/Moritz_body/data/body';
clc
a = dir('/Volumes/Moritz_body/data/body/**/*.avi');
fid = fopen('paths.txt','w');
com_bef = 'ffmpeg -i "';
com_mid = '" -c:a copy -preset ultrafast "';



for k = 1:length(a)
    path = strcat(a(k).folder,'/',a(k).name);
    com_af = strcat(a(k).folder,'/',a(k).name(1:end-4),'_r.avi"');
    if isempty(strfind(path,'/.')) &&isempty(strfind(path,'too'))
        fprintf(fid,[com_bef,path,com_mid,com_af,'\n']);
    end
        
end
fclose(fid);
%% 
path = '/Volumes/Moritz_body/data/body/18-05-03/2905l/1/2905l-s1.mp4';
data = VideoReader(path);


%%
filepaths = dir('/Volumes/Moritz_body/data/body/**/*.mp4');
all = cell(length(filepaths),2);

for f = 1:length(filepaths)
    if ~strcmp(filepaths(f).name(1),'.') 
        strcat(filepaths(f).folder,'/',filepaths(f).name)
        data = VideoReader(strcat(filepaths(f).folder,'/',filepaths(f).name));
        all{f,1} = data.FrameRate;
        all{f,2} = strcat(filepaths(f).folder,'/',filepaths(f).name);
    end
end

%%

p = '/Volumes/Moritz_body/data/body/17-10-16/5212r/215/5212r-s215_r.avi';

vid = VideoReader(p);

