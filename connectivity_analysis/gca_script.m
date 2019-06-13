%% GC ANALYSIS SCRIPT

% Meta variables
ROOT = "/Users/jannisborn/Desktop/HIFO/";
DATA_PATH = strcat(ROOT, "widefield_for_network/");
DISKS = ["Moritz_wide", "Moritz_img2"];
BRAIN_AREAS = ["PL", "ACC", "M2", "M1", "S1b", "S1mn", "S1bf", "S2", ...
    "Au", "ASA", "V2L", "RL", "V2A", "V1", "RS"];
   
NUM_SAMPLES = 200; % samples in widefield data per trial
SR        = 20;    % sample rate (Hz)
MOMAX     = 6;     % maximum model order for model order estimation (300ms delay)


PHASES = ["sensory", "response", "reward"];
TIMES = {[30, 66], [76, 126], [120, 160]};


%% Task specific extractions

MICE = ["1113rr", "1110r","1111lr", "5212r", "5627rr"];
TASKS = ["visual_task", "sensory_task", "naive_task"];
STIMS = ["V", "S", "N", "V+S"];
BEHS = ["H", "M", "FA", "CR", "EL"];
PHASES = ["sensory", "response"];
TIMES = {[30, 66], [76, 126]};

%MICE = ["1113rr"]; TASKS = ["visual_task"]; BEHS = ["H"];
%%
mouse_ind = 0;
try
    for mouse = MICE
        
        mouse_ind = mouse_ind + 1;
        disp(strcat("Starting with mouse ", num2str(mouse_ind), "/", ...
            num2str(length(MICE)), "."));
        mouse_path = strcat(ROOT, "gc_results/", mouse);
        mkdir(mouse_path);
        for task = TASKS
            if ~(mouse=="1113rr" && task == "visual_task")
            for stim = STIMS   
                phase_ind=0;
                for phase=PHASES
                    phase_ind=phase_ind+1;
                    disp(strcat("Now starting with mouse = ", mouse, ...
                        " task = ", task, " stim = ", stim, " phase = ", phase));
                    
                    brain_areas = BRAIN_AREAS;
                    % Assemble data matrices
                     data = zeros(length(brain_areas), NUM_SAMPLES, 0);
                     for disk = DISKS
                        file_path = strcat(DATA_PATH,"gca_",task,"_",mouse,"_disk_",disk);

                        % Load data from specific disk
                        try
                            d = data_loader_task_specific(file_path, task, stim, "A");
                            data = cat(3, data, d);
                        catch
                            warning(strcat("Path ", file_path, " not available"));
                        end

                    end
                    if size(data,3) == 0
                        disp("No data, skipping this computation");
                        continue
                    end

                    % Throw out missing ROIs
                    tmp = reshape(data, [size(data,1),size(data,2)*size(data,3)]);
                    missing_rois = find(~any(tmp,2));

                    if ~isempty(missing_rois)

                        data(missing_rois,:,:) = [];
                        brain_areas(missing_rois) = [];
                        disp(strcat('Excluded ROI numbers = ', num2str(missing_rois), ' due to missing data.'));
                    end

                    % Chop data to
                    data = data(:, TIMES{phase_ind}(1):TIMES{phase_ind}(2), :);
                    save_path = strcat(ROOT, "gc_results/", mouse,"/",...
                        task,"_stim_",stim, "_", phase, "_phase/");
                    mkdir(save_path);
                    size(data)
                    

                    results = gc_analysis(data, SR, MOMAX, save_path, brain_areas);
                    
                    if results.info.error
                        disp("Finding broken trial...")
                        % Verify spectral radius of data and exclude broken trials
                        trials=size(data,3);
                        counter=0;
                        data_copy = [];
                        for t=1:size(data,3)
                            data_copy = cat(3,data_copy, data(:,:,t));

                            pass = test_spectral_radius(data_copy);
                            if ~pass
                                data_copy(:,:,end) = [];
                                counter=counter+1;
                            end
                            %disp([t,counter]);
                        end
                        data = data_copy;
                        disp(strcat("Of ", num2str(trials), " trials, ", num2str(counter), ...
                            " were broken."));
                        results = gc_analysis(data, SR, MOMAX, save_path, brain_areas);
                    end
                end
            end
            end
            
            if task~="naive_task"
                for beh=BEHS
                    phase_ind=0;
                    for phase=PHASES
                        phase_ind=phase_ind+1;
                        disp(strcat("Now starting with mouse = ", mouse, ...
                            " task = ", task, " beh = ", beh, " phase = ", phase));

                        brain_areas = BRAIN_AREAS;
                        % Assemble data matrices
                         data = zeros(length(brain_areas), NUM_SAMPLES, 0);
                         for disk = DISKS
                            file_path = strcat(DATA_PATH,"gca_",task,"_",mouse,"_disk_",disk);

                            % Load data from specific disk
                            try
                                d = data_loader_task_specific(file_path, task, "A", beh);
                                data = cat(3, data, d);
                            catch
                                warning(strcat("Path ", file_path, "not available"));
                            end

                        end
                        if size(data,3) == 0
                            disp("No data, skipping this computation");
                            continue
                        end

                        % Throw out missing ROIs
                        tmp = reshape(data, [size(data,1),size(data,2)*size(data,3)]);
                        missing_rois = find(~any(tmp,2));

                        if ~isempty(missing_rois)

                            data(missing_rois,:,:) = [];
                            brain_areas(missing_rois) = [];
                            disp(strcat('Excluded ROI numbers = ', num2str(missing_rois), ' due to missing data.'));
                        end

                        % Chop data to
                        data = data(:, TIMES{phase_ind}(1):TIMES{phase_ind}(2), :);
                        save_path = strcat(ROOT, "gc_results/", mouse,"/",...
                            task,"_beh_",beh, "_", phase, "_phase/");
                        mkdir(save_path);
                        results = gc_analysis(data, SR, MOMAX, save_path, brain_areas);
                        
                        if results.info.error
                            disp("Finding broken trial...")
                            trials=size(data,3);
                            % Verify spectral radius of data and exclude broken trials
                            counter=0;
                            data_copy = [];
                            for t=1:size(data,3)
                                data_copy = cat(3,data_copy, data(:,:,t));

                                pass = test_spectral_radius(data_copy);
                                if ~pass
                                    data_copy(:,:,end) = [];
                                    counter=counter+1;
                                end
                                %disp([t,counter]);
                            end
                            data = data_copy;
                            disp(strcat("Of ", num2str(trials), " trials, ", num2str(counter), ...
                                " were broken."));
                            results = gc_analysis(data, SR, MOMAX, save_path, brain_areas);
                        end  
                                            
                    end
                end
            end
        end
    end
catch
    warning("Error in gca_script.m");
end


%% Fixing specific task-specific runs (those with non-stationary trials)
MICE = ["1113rr", "1110r", "1111lr", "5212r", "5627rr"];
TASKS = ["visual_task", "sensory_task", "naive_task"];
STIMS = ["V", "S", "N", "V+S"];
BEHS = ["H", "M", "FA", "CR", "EL"];
PHASES = ["sensory", "response"];
TIMES = {[30, 66], [76, 126]};


%%
for ind = 1:length(MICE)
    mouse =  MICE(ind);
    task = TASKS(ind);
    stim = STIMS(ind);
    phase = PHASES(ind);
    typ = TYPES(ind);
    
    if phase=="sensory"
        phase_ind = 1;
    elseif phase=="response"
        phase_ind = 2;
    end
    
    
    disp(strcat("Mouse ", mouse, ", task ", task, " stim, ", stim));
    mouse_path = strcat(ROOT, "gc_results/", mouse);
    mkdir(mouse_path);
        
    brain_areas = BRAIN_AREAS;
    % Assemble data matrices
    data = zeros(length(brain_areas), NUM_SAMPLES, 0);
    for disk = DISKS
        file_path = strcat(DATA_PATH,"gca_",task,"_",mouse,"_disk_",disk);

        % Load data from specific disk
        try
            if typ == "s"
                d = data_loader_task_specific(file_path, task, stim, "A");
                data = cat(3, data, d);
                s = "_stim_";
            elseif typ == "b"
                d = data_loader_task_specific(file_path, task, "A", stim);
                data = cat(3, data, d);
                s = "_beh_";
            end
        catch
            warning(strcat("Path ", file_path, " not available"));
        end

    end
    if size(data,3) == 0
        disp("No data, skipping this computation");
        continue
    end

    % Throw out missing ROIs
    tmp = reshape(data, [size(data,1),size(data,2)*size(data,3)]);
    missing_rois = find(~any(tmp,2));

    if ~isempty(missing_rois)

        data(missing_rois,:,:) = [];
        brain_areas(missing_rois) = [];
        disp(strcat('Excluded ROI numbers = ', num2str(missing_rois), ' due to missing data.'));
    end

    % Chop data to
    data = data(:, TIMES{phase_ind}(1):TIMES{phase_ind}(2), :);
    
    save_path = strcat(ROOT, "gc_results/", mouse,"/",...
        task,s,stim, "_", phase, "_phase/");
    mkdir(save_path);
    size(data)


    results = gc_analysis(data, SR, MOMAX, save_path, brain_areas);

    if results.info.error
        trials=size(data,3);
        disp("Finding broken trial...")
        % Verify spectral radius of data and exclude broken trials
        counter=0;
        data_copy = [];
        for t=1:size(data,3)
            data_copy = cat(3,data_copy, data(:,:,t));

            pass = test_spectral_radius(data_copy);
            if ~pass
                data_copy(:,:,end) = [];
                counter=counter+1;
            end
            %disp([t,counter]);
        end
        data = data_copy;
        disp(strcat("Of ", num2str(trials), " trials, ", num2str(counter), ...
            " were broken."));
        results = gc_analysis(data, SR, MOMAX, save_path, brain_areas);
    end
end


