
value = 1;
%app.meta.UI.            
close all;
img = read(vid.data,12);

drawnow
for roi = 1:length(vid.ROI.polygon)
    %polygon = reshape(vid.ROI.polygon{k},[len/2,2]);
    
    roi_points = vid.ROI.polygon{roi};
    for l = 1:length(roi_points)
         if roi == 1
            img(roi_points(l,1),roi_points(l,2),:) = [255 0 0]; % Mark it red
        elseif roi == 2
            img(roi_points(l,1),roi_points(l,2),:) = [0 255 0]; % Mark it green
        elseif roi == 3
            img(roi_points(l,1),roi_points(l,2),:) = [0 0 255]; % Mark it blue
        elseif roi ==4
            img(roi_points(l,1),roi_points(l,2),:) = [255 0 255]; % Mark it magenta
         end 
    end
end
    imshow(img)

%%
clc
close all;
ma = 1; % Pearson R motion method
trial = 3
figure
for roi = 1:size(motion.Results,1)
    vec = squeeze(motion.Results(roi,ma,motion.bounds(trial,1):motion.bounds(trial,2)));
    vec = (vec - min(vec)) / (max(vec) - min(vec)); % normalize to [0,1]
    
    
    if roi == 1 
        plot(1:length(vec),vec,'r'); hold on;
    elseif roi == 2
        plot(1:length(vec),vec+1.5,'g'); hold on;
    elseif roi == 3
        plot(1:length(vec),vec+3,'b'); hold on;
    elseif roi ==4
        plot(1:length(vec),vec+4.5,'m'); hold on;
    end 
    
    
end
y_uplim = 1+(roi-1)*1.5;
yticks(0:0.5:y_uplim); 
lab = [0:0.5:1];
for k = 1:roi-1
    lab = [lab,lab(1:3)];
end
lab
yticklabels(lab)
ylim([0 y_uplim]);
l1 = line([50 50], [0 y_uplim],'Color','black'); % how to delete line
l2 = line([100 100], [0 y_uplim],'Color','black');
delete(l1)


%%
ma = 'Pearson R';
for algo = 1:length(motion.algos)
    if strcmp(motion.algos{algo},ma)
end
                

    

