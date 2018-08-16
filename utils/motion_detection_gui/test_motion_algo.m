function [outputArg1,outputArg2] = test_motion_algo(path,trial,roi)
% This function verifies the functionality of the motion detection
% algorithms.
%
% Parameters:
% --------------
% PATH      {str} full path to folder with .mp4 video
% TRIAL     {int} the trial that should be verified
% ROI       {int} the ROI that is evaluated


% Retrieve data
vid = VideoReader(strcat(path,'/chunked/Trial_',num2str(trial),'.mp4'));
load(strcat(path,'/motion_analysis'));

strcat("We evaluate the ROI ", motion.ROIs{roi})

% lower bound is included, upper not anymore, thus no + 1
num_frames = motion.bounds(trial,2) - motion.bounds(trial,1);

% Prepare figure
figure
subplot(2,3,1);imshow(zeros(300)); title('Video');
subplot(2,3,2);imshow(zeros(300)); title('AbsDiff');
subplot(2,3,3);imshow(zeros(300)); title('Pearson R');
subplot(2,3,4);imshow(zeros(300)); title('Lucas Kanade');
subplot(2,3,5);imshow(zeros(300)); title('Horn Schunck');
titles = {'AbsDiff', 'Pearson R', 'Lucas Kanade', 'Horn Schunck'};

% Invert Pearson Results so that it makes sense visually
motion.Results(:,2,:) = (motion.Results(:,2,:) - 1) * (-1);

% To normalize the plots
norms = max(motion.Results(roi,:,:),[],3);

counter = 1;
for f = motion.bounds(trial,1):motion.bounds(trial,2)-1
    subplot(2,3,1);imshow(read(vid, counter)); title('Video');
    for k = 1:4
        subplot(2,3,k+1); imshow(motion.Results(roi,k,f)/norms(k) * ones(300));
        title(titles{k});
    end
    counter = counter + 1;
    drawnow
end

close all;

end