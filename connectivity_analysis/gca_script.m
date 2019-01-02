%% GC ANALYSIS SCRIPT

% Meta variables
ROOT = "/Users/jannisborn/Desktop/HIFO/";
DATA_PATH = strcat(ROOT, "gc_matrices/");
DISKS = ["F", "I"];
BRAIN_AREAS = ["V1 rostral", "V1 middle lateral", "V1 middle medial", ...
       "V1 caudal lateral", "V1 caudal medial", "V rostral", "V anterolateral", ...
       "V lateral", "V posteromedial", "V anterior", "Retrosplenial caudal", ...
       "Retrosplenial rostral", "Retrosplenial inter", "Aud", "S1 rostral", ...
       "S1 lateral", "S1 medial", "S1 caudal", "S trunk", "S hindlimb", ...
       "S forelimb", "S nose", "S mouth", "S supplementary", "M1 lateral", ...
       "M1 intermediate", "M1 medial", "M2 rostrolateral", "M2 rostromedial", ...
       "M2 intermiedate", "M2 caudal", "PFC prelimbic", "PFC ACC"];
   
   
SR        = 20;    % sample rate (Hz)
MOMAX     = 6;     % maximum model order for model order estimation (300ms delay)
% widefield data has 200 samples (10 sec) per trial. By default all samples
% are used (onset=1). Since 31 is stimulus onset, it may be benefical to try
% onset=31 to discard measurements related 
TIME_ONSET = 1;  


   
% Loop variables
% Data for the following mice is available: 5212r,1110r, 5627rr, 1113rr.
MOUSES = ["5212r"];
% Data for the following conditions is available: hit, correct_rejection,
%   miss, false_alarm, early_lick, multi_stim, no_stim, visual_stim,
%   sensory_stim.
CONDITIONS = ["hit"];
%% Loop over all mice and all conditions

try
    for condition = CONDITIONS
        for mouse = MOUSES

            
            disp(strcat("Now starting with mouse = ", mouse, " condition = ", condition));

            brain_areas = BRAIN_AREAS; % create local copy

            % Collect the data matrices
            data = zeros(33, 200, 0);
            for disk = DISKS
                file_path = strcat(DATA_PATH,"gca_",condition,"_",mouse,"_disk_",disk);

                % Load data from specific disk
                try
                    load(file_path);
                    data = cat(3, data, eval(strcat(condition, ".data")));
                catch
                    warning(strcat("Path ", file_path, "not available"));
                end
                
            end
            if size(data,3) == 0
                disp("Now data, skipping this computation");
                continue
            end

            % Throw out missing ROIs
            tmp = reshape(data, [size(data,1),size(data,2)*size(data,3)]);
            missing_rois = find(~any(tmp,2));
            data(missing_rois,:,:) = [];
            brain_areas(missing_rois) = [];
            disp(strcat('Excluded ROI numbers = ', num2str(missing_rois), ' due to missing data.'));

             % Remove sparsity bug
             % Should not be necessary anymore when data is recreated
             % properly.
            if ~(condition == "visual_task" || condition == "sensory_task")
                new_data = zeros(size(data));
                ind_old = 1;
                ind_new = 1;
                while ind_old < size(data,3)
                    roi_ind = 0;
                    slice = data(:,:,ind_old);
                    roi_ind_new = find(any(slice,2));
                    while roi_ind_new > roi_ind
                        roi_ind = roi_ind_new;
                        new_data(roi_ind,:,ind_new) = data(roi_ind, :, ind_old);
                        ind_old = ind_old+1;
                        if ind_old > size(data,3) continue; end
                        slice = data(:,:,ind_old);
                        roi_ind_new = find(any(slice,2));
                        
                    end
                    if isempty(roi_ind_new)
                        ind_old = ind_old + 1;
                    end
                    ind_new = ind_new + 1;
                    
                end
                new_data(:,:,ind_new+1:end) = [];
                data = new_data;
                disp("Done with removing sparsity bug");
            end
            
            % Throw out first X measurements if applicable
            data = data(:, TIME_ONSET:end, :);

            save_path = strcat(ROOT, "gc_results/", mouse,"_",condition,"/");
            mkdir(save_path);
            results = gc_analysis(data, SR, MOMAX, save_path, brain_areas);

        end

    end

catch
    warning("Error in gca_script.m");
end