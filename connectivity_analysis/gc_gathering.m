%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CONTEXT AND SOURCES:
% Class for data gathering for granger causality of Ca2+ data of multiple ROI 
%   of mice behaving under different task conditions. Data similar to:
%   [1] "Context-dependent cortical integration of visual and somatosensory 
%           stimuli in behaving mice" M.Buchholz, Y.Sych, F.Helmchen, A.Ayaz
% 
% FUNCTION:
% This is a preparatory of the gc_analysis class. 
% It gathers CA imaging data of all sessions of a given mouse that match
% a given criterion (e.g. visual task). It expects to receive 2 criteria, 
% retrieves the matching sessions/trials and bundles them into 2 data 
% matrices of shape 
%       #ROI x samples_per_trial x num_trial 
% that are saved within a struct and that can be used for gc analysis.
% 
%   Jannis Born, October 2018

classdef gc_gathering

properties (Access = public)
    % Declare class variables
    mouse
    group_crit
    gca % this is the struct that will be exported.
    
end

properties (Access = private)
   % Internal class variables
    
   allowed_mice = ["5627rr", "5212r", "1110r", "1111lr", "1113rr", "2905l", "2907ll"];
   allowed_group_crits = ["sensory_task", "visual_task", "naive_task", ...
       "sensory_stim", "visual_stim", "multi_stim", "no_stim", ...
       "hit", "miss", "false_alarm", "correct_rejection", "early_lick"];
   allowed_time_lag = [0, 1000]; 
   
   % Path variables
   file_paths = [ "F:\data\registered", "I:\data\registered"];
   delimiter = '\'; % / for Mac, \ for Windows
   meta_file_root = "H:\data\behavior\";
  
   % Data locations
   affine = load('C:\Users\Ayaz-Studi2\Desktop\Jannis\transformation_matrices');
   rois = load('C:\Users\Ayaz-Studi2\Desktop\Jannis\all_ROIs');
   brain_areas = ["V1_rostral", "V1_middle_lateral", "V1_middle_medial", ...
       "V1_caudal_lateral", "V1_caudal_medial", "V_rostral", "V_anterolateral", ...
       "V_lateral", "V_posteromedial", "V_anterior", "Retrosplenial_caudal", ...
       "Retrosplenial_rostral", "Retrosplenial_inter", "Aud", "S1_rostral", ...
       "S1_lateral", "S1_medial", "S1_caudal", "S_trunk", "S_hindlimb", ...
       "S_forelimb", "S_nose", "S_mouth", "S_supplementary", "M1_lateral", ...
       "M1_intermediate", "M1_medial", "M2_rostrolateral", "M2_rostromedial", ...
       "M2_intermiedate", "M2_caudal", "PFC_prelimbic", "PFC_ACC"];
   
   ca_img_size = [256, 256];
   samples_per_trial = 200;
   warper;
   num_trials;% Tracking number of trials of current session

end


methods (Access = public)

    % Constructor. 
    %
    % Class should be initialized by the following vars:
    % (1) mouse -  {<string>}, specifying name of the mouse. Choose 
    %           from {"5627rr", "5212r", "1110r", "1111lr", "1113rr"
    %           "2905l", "2907ll"}. Mouse ID for which analysis is
    %           performed.
    % (2) group_crit - {<string array>} of shape 1x2. Each element should
    %           be choosen from {"V", "S", "N"}. The data fed to the gc
    %           analysis is gathered from each session/trial fullfiling 
    %           the specified condition. The output of both GC analyses
    %           is then tested for significance.
    %           TODO: Enable further groupings such as 
    %               (A)     Stimulus type (V, S, N, VS).
    %               (B)     Response type (H, FA, CR, M).
    %               (C)     Performance level (High, low).
    % DEPRECATED
    % (3) max_time_lag - {<double>, <intX>}. Specifies the maximal time 
    %           lag (maximum model order) for the GC analysis. How much
    %           time is maximally allowed for the signal to travel from
    %           region A to B. Give in ms. Defaults to 300ms, corres-
    %           ponding to 15 frames (SR: 20 Hz). max_time_lag is only
    %           an upper bound, ideal model order selected via BIC
    % (4) padding - {<string>} from {"VALID", "SAME"}. [2] expects a
    %           matrix #ROI x frames_per_trial x #trials. In [1]
    %           frames_per_trial = 200 for all data, but #trials varies
    %           highly from session to session

    function obj = gc_gathering(varargin)

        % Error handling
        if (nargin < 2 || ~isstring(varargin{1}) || ~isstring(varargin{2})...
                || length(varargin{2})~= 2)
            error(['Please ensure the first arg is a STRING for the '...
                'mouse name and second arg is a STRING of length 2 '...
                'for the 2 groups to compare']);

        elseif nargin > 2
            warning("Third and all later args are discarded");
        end

        if any(contains(obj.allowed_mice, varargin{1}))
            obj.mouse = varargin{1};
        else
            error("Unknown mouse name given")
        end

        if sum(contains(obj.allowed_group_crits, varargin{2})) == 2 
            obj.group_crit = varargin{2};
        else
            error("Unknown or too many group criteria given")
        end

        
        %affine2d class to warp to standard atlas
        obj.warper = eval(strcat('obj.affine.transform.mouse',obj.mouse)); 
        
    end


    function assemble_data(obj)
    % This function assembles the data matrices. It searches for all
    % sessions of a given mouse and extracts the trials according to
    % the criteria defined in group_crit.
    
    % Allocate struct that will be exported
    obj.allocate_outputs();
    
    trial_sum = 10000; % just for data array allocation
    gc_ind = 0;
    
    
    for gc = obj.group_crit
        
        gc_ind = gc_ind + 1;
        data = zeros(length(obj.brain_areas), obj.samples_per_trial, trial_sum, 2);
        cti = 0; % cumulative_trial_index
        sessions = [];
        trials = [];
        
        % For each hard-disk, go through all subfolders
        for file_path = obj.file_paths
            
            disp(strcat("Searching now on disk ", file_path(1)));
            date_folders = obj.list_subfolders(file_path);

            % Open first level of folder (folder names are dates)
            for folder_ind = 1:length(date_folders)
                
                disp(strcat("Folder ", num2str(folder_ind), " out of ", ...
                    num2str(length(date_folders))));
                date_folder = date_folders(folder_ind);
                tmp = strsplit(date_folder,obj.delimiter);
                date_id = tmp(end);
                mouse_folders = obj.list_subfolders(date_folder);

                % Open 2nd level (folder names are mouse IDs)
                for mouse_folder = mouse_folders
                    tmp = strsplit(mouse_folder,obj.delimiter);
                    folder_name = tmp(end);

                    if strcmpi(folder_name, obj.mouse)
                        session_folders = obj.list_subfolders(mouse_folder);

                        % Open 3rd level (folder names are session IDs)
                        for session_folder = session_folders
                            tmp = strsplit(session_folder,obj.delimiter);
                            session_id = tmp(end);

                            % Check whether the session is of given type
                            % TODO: If other group criteria are allowed,
                            % this needs modification (optional execution
                            % only)
                            meta_path = strcat(obj.meta_file_root, date_id, ...
                                obj.delimiter, obj.mouse, obj.delimiter, ...
                                session_id, obj.delimiter);
                            try
                              metastats = readExperimentData(meta_path, ...
                                  strcat(obj.mouse,'-s',session_id,...
                                  '-exp.txt'));
                            catch
                                warning(strcat("Date ", date_id, ...
                                    " session ", session_id, "is ", ...
                                    "skipped, because metafile was ", ...
                                    "not found or threw an error"));

                            end
                            
                            session_type = obj.get_session_type(metastats);
                            
                            if strcmp(session_type, gc)
                                disp(['Now processing session ', ...
                                    num2str(session_id), ' recorded at ', ....
                                    num2str(date_id), '.']);
                                
                                t = obj.parse_data(session_folder);
                                data(:,:,cti+1 : cti+size(t,3)) = t;
                                cti = cti + size(t,3);
                                sessions = [sessions, session_id];
                                trials = [trials, obj.num_trials];
                            end
                        end
                    end
                end
            end
        end
        
        data(:,:,cti+1:end) = []; % remove unused array space.
        if gc_ind == 1
            obj.gca.data_1.data = data;
            obj.gca.data_1.sessions = sessions;
            obj.gca.data_1.trials = trials;
        elseif gc_ind == 2
            obj.gca.data_2.data = data;
            obj.gca.data_2.sessions = sessions;
            obj.gca.data_2.trials = trials;
        end
    end
    
    % descriptive name
    save_path = strcat(pwd, '/data/gca_',obj.group_crit(1), '_vs_', ...
        obj.group_crit(2), '_', obj.mouse);
    save(save_path,'obj.gca');
    
    end

end


methods (Access = private) % Internal methods

    function date_folders = list_subfolders(obj, path)
        % Receives a path to a folder and returns a cell array of sub-
        % folder-paths.
        sub_files = dir(path);
        dir_inds = [sub_files(:).isdir];
        date_folders = {sub_files(dir_inds).name}';
        date_folders(startsWith(date_folders,'.')) = [];
        date_folders = strcat(path, obj.delimiter, date_folders)';
    end
    
    function session_type = get_session_type(obj, stats)
        % Receives the stats of a session and returns whether the session
        % was visual ("V"), somatosensory ("S") or naive ("N") --> ?.
        
        if sum(strcmpi(stats.stim,"S") & strcmpi(stats.beh,"H")) > 0 
            session_type = "sensory_task";
        elseif sum(strcmpi(stats.stim,"V") & strcmpi(stats.beh,"H")) > 0 
            session_type = "visual_task";
        else
            warning(strcat("Session skipped, since type unclear "));
            session_type = "unknown";
        end
        
    end
    
    function allocate_outputs(obj)
        % Allocate and export a struct object
        obj.gca = struct();
        obj.gca.mouse = obj.mouse;
        obj.gca.rois = obj.brain_areas;
        obj.gca.groupings = obj.group_crit;
        obj.gca.data_1 = struct();
        obj.gca.data_2 = struct();
    end
    
    function data = parse_data(obj, path)
        % Receives the path to a given session and the metastats of the
        % session. Loops over all trials, reads the CA2+ matrix, warps the 
        % data to the standard atlas, applies the ROI mask for each
        % preserved ROI, averages and returns a matrix 'data' of shape:
        % #ROI x samples_per_trial x num_trials
        
        registration = load(strcat(path,'\registration'));
        
        % Allocations
        % Unfortunately data is double precision float (64bit), but still
        % in range of numerical imprecisions [0, 2]. Can't reliably convert
        % to 8bits (uint8) since upper bound unknown. 
        % Thus can't store data from all trials in memory since 
        %   256 x 256 x 200 x num_trials (450) are already ca. 45GB memory.
        
        %warped_data = zeros(obj.ca_img_size(1), obj.ca_img_size(2), ...
        %    obj.samples_per_trial, registration.info.trials_obj);
        data = zeros(length(obj.brain_areas), obj.samples_per_trial, ...
            registration.info.trials_obj);
        cutter = imref2d(obj.ca_img_size);
        
        % Tracking number of trials of current session
        obj.num_trials = registration.info.trials_obj;

        % Load imaging data trial per trial
        tic
        for trial = 1:registration.info.trials_obj
            if mod(trial, 20) == 0
                disp(['Currently processing trial ', num2str(trial),...
                    '/', num2str(registration.info.trials_obj)]);
                toc
                tic;
            end           
            trial_data = load(strcat(path,'\dFF_t',num2str(trial)));
            %warped_data(:,:,:,trial) = imwarp(trial_data.dFF, obj.warper, ...
            %    'OutputView', cutter);
            warped_data = imwarp(trial_data.dFF, obj.warper,'OutputView', cutter);
            warped_data = reshape(warped_data, [size(warped_data,1)*...
                size(warped_data,2),size(warped_data,3)]);
            
            brain_area_ind = 0;
            % Save average response of all rois for all frames of all trials.
            for brain_area = obj.brain_areas

                brain_area_ind = brain_area_ind+1;
                area_mask = eval(strcat('obj.rois.ROIs.',brain_area,'.maskCircle'));
                flat_mask = reshape(area_mask, [size(area_mask,1)*...
                    size(area_mask,2),1]);

                data(brain_area_ind,:,trial) = squeeze(mean(warped_data(flat_mask,:),1));
            end

        end
            
        %w = reshape(warped_data, [size(warped_data,1)*size(warped_data,2), ...
        %    size(warped_data,3), size(warped_data,4)]);
        %size(warped_data)
        
        % Allocations
        %data = zeros(length(obj.brain_areas), obj.samples_per_trial, ...
        %    registration.info.trials_obj);
%         brain_area_ind = 0;
%         % Save average response of all rois for all frames of all trials.
%         for brain_area = obj.brain_areas
%             
%             brain_area_ind = brain_area_ind+1;
%             area_mask = eval(['obj.rois.ROIs.',brain_area,',maskCircle']);
%             flat_mask = reshape(area_mask, [size(area_mask,1)*...
%                 size(area_mask,2),1]);
%             
%             data(brain_area_ind,:,:) = squeeze(mean(w(flat_mask,:,:),1));
%         end
    end
end
    
end
