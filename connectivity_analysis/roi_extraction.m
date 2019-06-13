%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CONTEXT AND SOURCES:
% Class for data gathering for granger causality of Ca2+ data of multiple ROI 
%   of mice behaving under different task conditions. 
%   Data similar to:
%   [1] "Context-dependent cortical integration of visual and somatosensory 
%           stimuli in behaving mice" M.Buchholz, Y.Sych, F.Helmchen, A.Ayaz
% 
% FUNCTION:
% This is a preparatory of the gc_analysis class. 
% INPUTS:
%   (1) mouse -  {<string>}, specifying name of the mouse. Choose from
%       {"5627rr", "5212r", "1110r", "1111lr", "1113rr", "2905l", "2907ll"}.
%       Mouse ID for which analysis is performed.
%   (2) disk - {<string>, <char>}, specifying the disk to search on. 
%       Choose from {"I", "F"}
%
% SAVES:
%   For each of the following groups (see below), one struct is saved in
%       H:\Jannis\granger_causality_data\
%
%   Each file has shape #ROI x samples_per_trial x num_trial 
%
%   NOTE:
%      Sessions where the mouse was unreactive (always or never licked) 
%           are discarded.
% 
% groups = {sensory_task, visual_task, naive_task, sensory_stim, visual_stim, 
%   multi_stim, no_stim, hit, miss, false_alarm, correct_rejection, early_lick}
%
%   Jannis Born, October 2018

classdef roi_extraction

properties (Access = public)
    % Declare class variables
    mouse
    disk
    sensory_task; visual_task; naive_task; sensory_stim; visual_stim; 
    multi_stim; no_stim; hit; miss; false_alarm; correct_rejection;
    early_lick;
end

properties (Access = private)
   % Internal class variables
   allowed_mice = ["5627rr", "5212r", "1110r", "1111lr", "1113rr", "2905l", "2907ll", "2906r"];
   allowed_disks = ["Moritz_wide", "Moritz_img2"];

   % Path variables
   data_root_pre = "/Volumes/"
   data_root_post = "/data/registered";
   delimiter = '/'; % / for Mac, \ for Windows
   meta_file_root = "/Volumes/Moritz_beh/data/behavior/";
   %save_path = "H:\Jannis\granger_causality_data\";
   save_path = "/Volumes/Moritz_beh/Jannis/widefield_for_network/";
  
   % Data locations
   affine = load('/Volumes/Moritz_beh/Jannis/transformation_matrices');

   rois = load('/Volumes/Moritz_beh/Jannis/new_rois');
   brain_areas = ["PL","ACC","M2","M1","S1b","S1mn","S1bf","S2","Au", ...
       "ASA","V2L","RL","V2A","V1","RS"];
%    rois = load('/Volumes/Moritz_beh/Jannis/all_ROIs');
%    brain_areas = ["V1_rostral", "V1_middle_lateral", "V1_middle_medial", ...
%        "V1_caudal_lateral", "V1_caudal_medial", "V_rostral", "V_anterolateral", ...
%        "V_lateral", "V_posteromedial", "V_anterior", "Retrosplenial_caudal", ...
%        "Retrosplenial_rostral", "Retrosplenial_inter", "Aud", "S1_rostral", ...
%        "S1_lateral", "S1_medial", "S1_caudal", "S_trunk", "S_hindlimb", ...
%        "S_forelimb", "S_nose", "S_mouth", "S_supplementary", "M1_lateral", ...
%        "M1_intermediate", "M1_medial", "M2_rostrolateral", "M2_rostromedial", ...
%        "M2_intermiedate", "M2_caudal", "PFC_prelimbic", "PFC_ACC"];
   
   cutter = imref2d([256, 256]); % image size = 256 x 256
   samples_per_trial = 200;
   warper;
   num_trials;% Tracking number of trials of current session
   trial_sum = 10000; % just for data array allocation
   cti; % order of ctis is like in 'groups'.
   
   % data structs
   groups = 12; % count of all following groups.
end


methods (Access = public)

    % Constructor. 
    %
    % Class should be initialized by the following vars:
    % (1) mouse -  {<string>}, specifying name of the mouse. Choose 
    %           from {"5627rr", "5212r", "1110r", "1111lr", "1113rr"
    %           "2905l", "2907ll"}. Mouse ID for which analysis is
    %           performed.
    % DEPRECATED
    % (2) group_crit - {<string array>} of shape 1x2. Each element should
    %           be choosen from {"V", "S", "N"}. The data fed to the gc
    %           analysis is gathered from each session/trial fullfiling 
    %           the specified condition. The output of both GC analyses
    %           is then tested for significance.
    %            further groupings such as 
    %               (A)     Stimulus type (V, S, N, VS).
    %               (B)     Response type (H, FA, CR, M).
    %               (C)     Performance level (High, low).
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
        if nargin < 2 || ~isstring(varargin{1}) || ~isstring(varargin{2})
            error(['Please ensure the first arg is a STRING for the '...
                'mouse name and the second one for the disk name.']);
        elseif nargin > 2
            warning("Second and all later args are discarded");
        end

        if any(contains(obj.allowed_mice, varargin{1}))
            obj.mouse = varargin{1};
        else
            error(strcat("Unknown mouse name given ", varargin{1}, ...
                "check help for details."));
        end
        
        %if ~strcmpi(varargin{2},"I") && ~strcmpi(varargin{2}, "F")
        if any(contains(obj.allowed_disks, varargin{2}))
            obj.disk = varargin{2};
        else
             error(strcat("Unknown disk name given (", varargin{2}, ...
                ") see help for details.")); 
        end

        %affine2d class to warp to standard atlas
        %obj.warper = eval(strcat('obj.affine.transform.mouse',obj.mouse)); 
        
    end


    function assemble_data(obj)
    % This function assembles the data matrices. It searches for all
    % sessions of a given mouse and extracts the trials according to
    % the criteria defined in group_crit.
    
    % Allocate structs that will be exported
    obj = obj.allocate_outputs();
    % cumulative_trial_indices. Tracks for each group how many columns were
    % written.
    obj.cti = zeros(obj.groups,1); 

    disp(strcat("Searching on disk ", obj.disk, " for mouse ", obj.mouse));
    file_path = strcat(obj.data_root_pre, obj.disk, obj.data_root_post);
    date_folders = obj.list_subfolders(file_path);

    % Open first level of folder (folder names are dates)
    for folder_ind = 1:length(date_folders)

        %disp(strcat("Folder ", num2str(folder_ind), " out of ", ...
        %    num2str(length(date_folders))));
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
                    % Check whether folder name is valid session_id
                    [num, status] = str2num(session_id);
                    error = ~status;

                    % Check whether the session is of given type
                    if ~error
                        meta_path = strcat(obj.meta_file_root, date_id, ...
                            obj.delimiter, obj.mouse, obj.delimiter, ...
                            session_id, obj.delimiter);
                    try
                        metastats = readExperimentData(meta_path, ...
                          strcat(obj.mouse,'-s',session_id,...
                          '-exp.txt'));
                    catch
                        warning(strcat("Date ", date_id, ...
                            " session ", num2str(session_id), "is ", ...
                            "skipped, because metafile was ", ...
                            "not found or threw an error"));
                        error = 1;
                    end

                    %session_type = obj.get_session_type(metastats, date_id, session_id);
                    session_type = metastats.task;

                    %if ~strcmpi(session_type,"unknown") && ~error
                    if ~isempty(session_type) && ~error
                        disp(['Now processing session ',num2str(session_id),...
                            ' recorded at ', num2str(date_id), ' with ', ...
                            num2str(length(metastats.beh)), 'trials.']);
                        obj = obj.parse_data(session_folder, session_type, metastats);
                    end
                    end
                end
            end
        end
    end
    
    obj.save_all()  
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
    
    function session_type = get_session_type(obj, stats, date_id, session_id)
        % Receives the stats of a session and returns whether the session
        % was visual ("V"), somatosensory ("S") or naive ("N") --> ?.
        
        if sum(strcmpi(stats.stim,"S") & strcmpi(stats.beh,"H")) > 0 
            session_type = "sensory_task";
        elseif sum(strcmpi(stats.stim,"V") & strcmpi(stats.beh,"H")) > 0 
            session_type = "visual_task";
        else
            warning(strcat("Date ",date_id," session ",num2str(session_id), ...
                " is skipped, since type unclear "));
            session_type = "unknown";
        end
        
    end
    
    function ind = get_index(obj,strct)
        % Receive length of data struct and trial no and returns ind 
        % to write next data
        if length(fieldnames(strct.data)) == 0
            ind = 1;
        else
            ind = length(strct.data) + 1;
        end
%         if length(fieldnames(strct.data)) == 0
%             ind = length(fieldnames(strct.data)) + 1;
%         else
%             ind = l+1;
%         end
            
            
    end
    
    function obj = allocate_outputs(obj)
        % Allocate and export a struct object
        
        spt = obj.samples_per_trial;
        ts = obj.trial_sum;
        
    
        obj.sensory_task = struct();
        obj.sensory_task.mouse = obj.mouse;
        obj.sensory_task.rois = obj.brain_areas;
        obj.sensory_task.data = struct();
        
%         obj.sensory_task.data = zeros(length(obj.brain_areas), spt , ts);
%         obj.sensory_task.beh = cell(ts,1);
%         obj.sensory_task.stim = cell(ts,1);
%         obj.sensory_task.session = cell(ts,1);
        
        obj.visual_task = struct();
        obj.visual_task.mouse = obj.mouse;
        obj.visual_task.rois = obj.brain_areas;
        obj.visual_task.data = struct();
        
%         obj.visual_task.data = zeros(length(obj.brain_areas), spt , ts);
%         obj.visual_task.beh = cell(ts,1);
%         obj.visual_task.stim = cell(ts,1);
%         obj.visual_task.session = cell(ts,1);

        obj.naive_task = struct();
        obj.naive_task.mouse = obj.mouse;
        obj.naive_task.rois = obj.brain_areas;
        obj.naive_task.data = struct();
        
%         obj.naive_task.data = zeros(length(obj.brain_areas), spt , ts);
%         obj.naive_task.beh = cell(ts,1);
%         obj.naive_task.stim = cell(ts,1);
%         obj.naive_task.session = cell(ts,1);      

        obj.sensory_stim = struct();
        obj.sensory_stim.mouse = obj.mouse;
        obj.sensory_stim.rois = obj.brain_areas;
        obj.sensory_stim.data = struct();
        
%         obj.sensory_stim.data = zeros(length(obj.brain_areas), spt , ts);   
%         obj.sensory_stim.beh = cell(ts,1);
%         obj.sensory_stim.stim = cell(ts,1);
%         obj.sensory_stim.session = cell(ts,1);
        
        obj.visual_stim = struct();
        obj.visual_stim.mouse = obj.mouse;
        obj.visual_stim.rois = obj.brain_areas;
        obj.visual_stim.data = struct();

%         obj.visual_stim.data = zeros(length(obj.brain_areas), spt , ts);   
%         obj.visual_stim.beh = cell(ts,1);
%         obj.visual_stim.stim = cell(ts,1);
%         obj.visual_stim.session = cell(ts,1);
        
        obj.multi_stim = struct();
        obj.multi_stim.mouse = obj.mouse;
        obj.multi_stim.rois = obj.brain_areas;
        obj.multi_stim.data = struct();

%         obj.multi_stim.data = zeros(length(obj.brain_areas), spt , ts);  
%         obj.multi_stim.beh = cell(ts,1);
%         obj.multi_stim.stim = cell(ts,1);
%         obj.multi_stim.session = cell(ts,1); 
        
        obj.no_stim = struct();
        obj.no_stim.mouse = obj.mouse;
        obj.no_stim.rois = obj.brain_areas;
        obj.no_stim.data = struct();

%         obj.no_stim.data = zeros(length(obj.brain_areas), spt , ts);  
%         obj.no_stim.beh = cell(ts,1);
%         obj.no_stim.stim = cell(ts,1);
%         obj.no_stim.session = cell(ts,1);        
        
        obj.early_lick = struct();
        obj.early_lick.mouse = obj.mouse;
        obj.early_lick.rois = obj.brain_areas;
        obj.early_lick.data = struct();
        
%         obj.early_lick.data = zeros(length(obj.brain_areas), spt , ts);   
%         obj.early_lick.beh = cell(ts,1);
%         obj.early_lick.stim = cell(ts,1);
%         obj.early_lick.session = cell(ts,1);
        
        obj.correct_rejection = struct();
        obj.correct_rejection.mouse = obj.mouse;
        obj.correct_rejection.rois = obj.brain_areas;
        obj.correct_rejection.data = struct();
%         obj.correct_rejection.data = zeros(length(obj.brain_areas), spt , ts);   
%         obj.correct_rejection.beh = cell(ts,1);
%         obj.correct_rejection.stim = cell(ts,1);
%         obj.correct_rejection.session = cell(ts,1);
        
        obj.false_alarm = struct();
        obj.false_alarm.mouse = obj.mouse;
        obj.false_alarm.rois = obj.brain_areas;
        obj.false_alarm.data = struct();
%         obj.false_alarm.data = zeros(length(obj.brain_areas), spt , ts);   
%         obj.false_alarm.beh = cell(ts,1);
%         obj.false_alarm.stim = cell(ts,1);
%         obj.false_alarm.session = cell(ts,1);
        
        obj.miss = struct();
        obj.miss.mouse = obj.mouse;
        obj.miss.rois = obj.brain_areas;
        obj.miss.data = struct();

%         obj.miss.data = zeros(length(obj.brain_areas), spt , ts);   
%         obj.miss.beh = cell(ts,1);
%         obj.miss.stim = cell(ts,1);
%         obj.miss.session = cell(ts,1);
        
        obj.hit = struct();
        obj.hit.mouse = obj.mouse;
        obj.hit.rois = obj.brain_areas;
        obj.hit.data = struct();
%         obj.hit.data = zeros(length(obj.brain_areas), spt , ts);   
%         obj.hit.beh = cell(ts,1);
%         obj.hit.stim = cell(ts,1);
%         obj.hit.session = cell(ts,1);
%         
    end
    
    function save_all(obj)
    
    % remove unused array space.
%     obj.cti(obj.cti==0) = 1; % safety in case one group was not found.
%     obj.visual_task.data(:,:,obj.cti(1):end) = []; 
%     obj.sensory_task.data(:,:,obj.cti(2):end) = []; 
%     obj.naive_task.data(:,:,obj.cti(3):end) = []; 
% 
%     obj.visual_stim.data(:,:,obj.cti(4):end) = []; 
%     obj.sensory_stim.data(:,:,obj.cti(5):end) = []; 
%     obj.multi_stim.data(:,:,obj.cti(6):end) = []; 
%     obj.no_stim.data(:,:,obj.cti(7):end) = []; 
%     
%     obj.hit.data(:,:,obj.cti(8):end) = []; 
%     obj.miss.data(:,:,obj.cti(9):end) = []; 
%     obj.false_alarm.data(:,:,obj.cti(10):end) = []; 
%     obj.correct_rejection.data(:,:,obj.cti(11):end) = []; 
%     obj.early_lick.data(:,:,obj.cti(12):end) = []; 
%     
%     obj.visual_task.beh(~cellfun('isempty',obj.visual_task.beh));
%     obj.visual_task.stim(~cellfun('isempty',obj.visual_task.stim));
%     obj.visual_task.session(~cellfun('isempty',obj.visual_task.session));
% 
%     obj.sensory_task.beh(~cellfun('isempty',obj.sensory_task.beh));
%     obj.sensory_task.stim(~cellfun('isempty',obj.sensory_task.stim));
%     obj.sensory_task.session(~cellfun('isempty',obj.sensory_task.session));
% 
%     obj.naive_task.beh(~cellfun('isempty',obj.naive_task.beh));
%     obj.naive_task.stim(~cellfun('isempty',obj.naive_task.stim));
%     obj.naive_task.session(~cellfun('isempty',obj.naive_task.session));
% 
%     obj.visual_stim.beh(~cellfun('isempty',obj.visual_stim.beh))
%     obj.visual_stim.stim(~cellfun('isempty',obj.visual_stim.stim))
%     obj.visual_stim.session(~cellfun('isempty',obj.visual_stim.session))
% 
%     obj.sensory_stim.beh(~cellfun('isempty',obj.sensory_stim.beh))
%     obj.sensory_stim.stim(~cellfun('isempty',obj.sensory_stim.stim))
%     obj.sensory_stim.session(~cellfun('isempty',obj.sensory_stim.session))
% 
%     obj.multi_stim.beh(~cellfun('isempty',obj.multi_stim.beh))
%     obj.multi_stim.stim(~cellfun('isempty',obj.multi_stim.stim))
%     obj.multi_stim.session(~cellfun('isempty',obj.multi_stim.session))
% 
%     obj.no_stim.beh(~cellfun('isempty',obj.no_stim.beh))
%     obj.no_stim.stim(~cellfun('isempty',obj.no_stim.stim))
%     obj.no_stim.session(~cellfun('isempty',obj.no_stim.session))
% 
%     obj.hit.beh(~cellfun('isempty',obj.hit.beh))
%     obj.hit.stim(~cellfun('isempty',obj.hit.stim))
%     obj.hit.session(~cellfun('isempty',obj.hit.session))
% 
%     obj.miss.beh(~cellfun('isempty',obj.miss.beh))
%     obj.miss.stim(~cellfun('isempty',obj.miss.stim))
%     obj.miss.session(~cellfun('isempty',obj.miss.session))
% 
%     obj.false_alarm.beh(~cellfun('isempty',obj.false_alarm.beh))
%     obj.false_alarm.stim(~cellfun('isempty',obj.false_alarm.stim))
%     obj.false_alarm.session(~cellfun('isempty',obj.false_alarm.session))
% 
%     obj.correct_rejection.beh(~cellfun('isempty',obj.correct_rejection.beh))
%     obj.correct_rejection.stim(~cellfun('isempty',obj.correct_rejection.stim))
%     obj.correct_rejection.session(~cellfun('isempty',obj.correct_rejection.session))
% 
%     obj.early_lick.beh(~cellfun('isempty',obj.early_lick.beh))
%     obj.early_lick.stim(~cellfun('isempty',obj.early_lick.stim))
%     obj.early_lick.session(~cellfun('isempty',obj.early_lick.session))
    
    
    % save all variables
    visual_task = obj.visual_task;
    sensory_task = obj.sensory_task;
    naive_task = obj.naive_task;

    visual_stim = obj.visual_stim;
    sensory_stim = obj.sensory_stim;
    multi_stim = obj.multi_stim;
    no_stim = obj.no_stim;
    
    hit = obj.hit;
    miss = obj.miss;
    false_alarm = obj.false_alarm;
    correct_rejection = obj.correct_rejection;
    early_lick = obj.early_lick;
    
    save(strcat(obj.save_path,'gca_visual_task_', obj.mouse, '_disk_', obj.disk),'visual_task','-v7.3');
    save(strcat(obj.save_path,'gca_sensory_task_', obj.mouse, '_disk_', obj.disk),'sensory_task','-v7.3');
    save(strcat(obj.save_path,'gca_naive_task_', obj.mouse, '_disk_', obj.disk),'naive_task','-v7.3');
    save(strcat(obj.save_path,'gca_visual_stim_', obj.mouse, '_disk_', obj.disk),'visual_stim','-v7.3');
    save(strcat(obj.save_path,'gca_sensory_stim_', obj.mouse, '_disk_', obj.disk),'sensory_stim','-v7.3');
    save(strcat(obj.save_path,'gca_multi_stim_', obj.mouse, '_disk_', obj.disk),'multi_stim','-v7.3');
    save(strcat(obj.save_path,'gca_no_stim_', obj.mouse, '_disk_', obj.disk),'no_stim','-v7.3');
    save(strcat(obj.save_path,'gca_hit_', obj.mouse, '_disk_', obj.disk),'hit','-v7.3');
    save(strcat(obj.save_path,'gca_miss_', obj.mouse, '_disk_', obj.disk),'miss','-v7.3');
    save(strcat(obj.save_path,'gca_false_alarm_', obj.mouse, '_disk_', obj.disk),'false_alarm','-v7.3');
    save(strcat(obj.save_path,'gca_correct_rejection_', obj.mouse, '_disk_', obj.disk),'correct_rejection','-v7.3');
    save(strcat(obj.save_path,'gca_early_lick_', obj.mouse, '_disk_', obj.disk),'early_lick','-v7.3');
    end
    
    function obj = parse_data(obj, path, session_type, metastats)
        % Receives the path to a given session and the metastats of the
        % session. Loops over all trials, reads the CA2+ matrix, warps the 
        % data to the standard atlas, applies the ROI mask for each
        % preserved ROI, averages and returns a matrix 'data' of shape:
        % #ROI x samples_per_trial x num_trials
        registration = load(strcat(path,obj.delimiter,'registration'));
   
        % Allocations
        % Unfortunately data is double precision float (64bit), but still
        % in range of numerical imprecisions [0, 2]. Can't reliably convert
        % to 8bits (uint8) since upper bound unknown. 
        % Thus can't store data from all trials in memory since 
        %   256 x 256 x 200 x num_trials (450) are already ca. 45GB memory.
        
%         if strcmpi(session_type, "visual_task")
%             session_data = obj.visual_task.data; session_ind = obj.cti(1);
%             session_beh = obj.visual_task.beh;
%             session_stim = obj.visual_task.stim;
%             session_session = obj.visual_task.session;
%         elseif strcmpi(session_type, "sensory_task")
%             session_data = obj.sensory_task.data; session_ind = obj.cti(2);
%             session_beh = obj.sensory_task.beh;
%             session_stim = obj.sensory_task.stim;
%             session_session = obj.sensory_task.session;
%         elseif strcmpi(session_type, "naive_task")
%             session_data = obj.naive_task.data; session_ind = obj.cti(3);
%             session_beh = obj.naive_task.beh;
%             session_stim = obj.naive_task.stim;
%             session_session = obj.naive_task.session;
%         end
        % Load imaging data trial per trial
        if ~(registration.info.trials_obj==length(metastats.beh))
            warning(strcat("Found CA data for ", num2str(registration.info.trials_obj), ...
                " trials, but behavioral data for ", num2str(length(metastats.beh)), ...
                ". Will neglect overhead trials."));
        end
        
         tic
         for trial = 1:min(registration.info.trials_obj,length(metastats.beh))
%         for trial = 1:length(metastats.beh)
            
            error = 0;
            if mod(trial, 1000) == 0
                toc
                disp(['Currently processing trial ', num2str(trial),...
                    '/', num2str(length(metastats.beh))]);
                tic
            end
            try 
                trial_data = matfile(strcat(path,obj.delimiter,'dFF_t',num2str(trial)));

                data = trial_data.dFF;
                %if any(isnan(data))
                if any(isnan(data(:)))
                    error = 1;
                    warning(strcat("Trial ", num2str(trial), " was ", ...
                       "skipped since it contains at least one NaN."));
                end
                
            catch
                strcat(path,obj.delimiter,'dFF_t',num2str(trial));
                warning(strcat("Trial ", num2str(trial), " was not found."));
                error = 1;
            end
            if ~error
                %warped_data = imwarp(trial_data.dFF, obj.warper,'OutputView', obj.cutter);
                %warped_data = trial_data.dFF;
                %warped_data = reshape(warped_data, [size(warped_data,1)*...
                %    size(warped_data,2),size(warped_data,3)]);
                
                data = reshape(data, [size(data,1)*...
                    size(data,2), size(data,3)]);
                
                
%                 session_ind = session_ind + 1;
%                 session_stim{session_ind} = metastats.stim{trial};
%                 session_beh{session_ind} = metastats.beh{trial};
                
                
%                 stim = metastats.stim{trial};
%                 beh = metastats.beh{trial};
%                 tmp = strsplit(path,obj.delimiter);
%                 session_id = tmp(end);
%                 session_session{session_ind} = session_id;
                
                % FIRST: Write information that holds for the entire trial
                %   Start with task data.
                
                if session_type == "visual_task"
                    ind = obj.get_index(obj.visual_task);
                    obj.visual_task.data(ind).wf = zeros(length(obj.rois.rois),...
                        obj.samples_per_trial);
                    obj.visual_task.data(ind).beh = metastats.beh{trial};
                    obj.visual_task.data(ind).stim= metastats.stim{trial};
                elseif session_type == "sensory_task"
                    ind = obj.get_index(obj.sensory_task);
                    obj.sensory_task.data(ind).wf = zeros(length(obj.rois.rois),...
                        obj.samples_per_trial);
                    obj.sensory_task.data(ind).beh = metastats.beh{trial};
                    obj.sensory_task.data(ind).stim= metastats.stim{trial};
                elseif session_type == "naive_task"
                    ind = obj.get_index(obj.naive_task);
                    obj.naive_task.data(ind).wf = zeros(length(obj.rois.rois),...
                        obj.samples_per_trial);
                    obj.naive_task.data(ind).beh = metastats.beh{trial};
                    obj.naive_task.data(ind).stim= metastats.stim{trial};
                else
                    warning(strcat("Unknown task type", session_type))
                end
                
                % Now write stim.data
                if metastats.stim{trial} == "V"
                    ind = obj.get_index(obj.visual_stim);
                    obj.visual_stim.data(ind).wf = zeros(length(obj.rois.rois),...
                        obj.samples_per_trial);
                    obj.visual_stim.data(ind).beh = metastats.beh{trial};
                    obj.visual_stim.data(ind).session = session_type; 
                elseif metastats.stim{trial} == "S"
                    ind = obj.get_index(obj.sensory_stim);
                    obj.sensory_stim.data(ind).wf = zeros(length(obj.rois.rois),...
                        obj.samples_per_trial);
                    obj.sensory_stim.data(ind).beh = metastats.beh{trial};
                    obj.sensory_stim.data(ind).session = session_type; 
                elseif metastats.stim{trial} == "V+S"
                    ind = obj.get_index(obj.multi_stim);
                    obj.multi_stim.data(ind).wf = zeros(length(obj.rois.rois),...
                        obj.samples_per_trial);
                    obj.multi_stim.data(ind).beh = metastats.beh{trial};
                    obj.multi_stim.data(ind).session = session_type;
                elseif metastats.stim{trial} == "N"
                    ind = obj.get_index(obj.no_stim);
                    obj.no_stim.data(ind).wf = zeros(length(obj.rois.rois),...
                        obj.samples_per_trial);
                    obj.no_stim.data(ind).beh = metastats.beh{trial};
                    obj.no_stim.data(ind).session = session_type;
                else
                	warning(strcat("Unknown stimulus ", metastats.stim{trial}))
                end 
                
                % Now write beh data
                if metastats.beh{trial} == "H"
                    ind = obj.get_index(obj.hit);
                    obj.hit.data(ind).wf = zeros(length(obj.rois.rois),...
                        obj.samples_per_trial);
                    obj.hit.data(ind).stim = metastats.stim{trial};
                    obj.hit.data(ind).session = session_type; 
                elseif metastats.beh{trial} == "M"
                    ind = obj.get_index(obj.miss);
                    obj.miss.data(ind).wf = zeros(length(obj.rois.rois),...
                        obj.samples_per_trial);
                    obj.miss.data(ind).stim = metastats.stim{trial};
                    obj.miss.data(ind).session = session_type; 
                elseif metastats.beh{trial} == "FA"
                    ind = obj.get_index(obj.false_alarm);
                    obj.false_alarm.data(ind).wf = zeros(length(obj.rois.rois),...
                        obj.samples_per_trial);
                    obj.false_alarm.data(ind).stim = metastats.stim{trial};
                    obj.false_alarm.data(ind).session = session_type;
                elseif metastats.beh{trial} == "CR"
                    ind = obj.get_index(obj.correct_rejection);
                    obj.correct_rejection.data(ind).wf = zeros(length(obj.rois.rois),...
                        obj.samples_per_trial);
                    obj.correct_rejection.data(ind).stim = metastats.stim{trial};
                    obj.correct_rejection.data(ind).session = session_type;
                elseif metastats.beh{trial} == "EL"
                    ind = obj.get_index(obj.early_lick);
                    obj.early_lick.data(ind).wf = zeros(length(obj.rois.rois),...
                        obj.samples_per_trial);
                    obj.early_lick.data(ind).stim = metastats.stim{trial};
                    obj.early_lick.data(ind).session = session_type;
                elseif metastats.beh{trial} == "error"
                    tmp=42;
                else
                	warning(strcat("Unknown behavior ", metastats.beh{trial}))
                end 
                % NOW
                % Loop over ROIs and save values.
                    
                % Save average response of all rois for all frames of all trials.
                for roi_ind = 1:length(obj.brain_areas)
                    %area_mask = eval(strcat('obj.rois.ROIs.',obj.brain_areas(roi_ind),'.maskCircle'));
                    area_mask = obj.rois.rois(roi_ind).mask;
                    flat_mask = reshape(area_mask, [size(area_mask,1)*...
                         size(area_mask,2),1]);

    
                    roi_value = squeeze(mean(data(flat_mask,:),1));

                    % Now write data to the right arrays. Start with session
                    %session_data(roi_ind,:,session_ind) = roi_value;
                    if session_type == "visual_task"
                        ind = length(obj.visual_task.data);
                        obj.visual_task.data(ind).wf(roi_ind,:) = roi_value;
                    elseif session_type == "sensory_task"
                        ind = length(obj.sensory_task.data);
                        obj.sensory_task.data(ind).wf(roi_ind,:) = roi_value;
                    elseif session_type == "naive_task"
                        ind = length(obj.naive_task.data);
                        obj.naive_task.data(ind).wf(roi_ind,:) = roi_value;
                    end
                    
                    % Now write stimulus data
                    if metastats.stim{trial} == "V"
                        ind = length(obj.visual_stim.data);
                        obj.visual_stim.data(ind).wf(roi_ind,:) = roi_value;
                    elseif metastats.stim{trial} == "S"
                        ind = length(obj.sensory_stim.data);
                        obj.sensory_stim.data(ind).wf(roi_ind,:) = roi_value;
                    elseif metastats.stim{trial} == "V+S"
                        ind = length(obj.multi_stim.data);
                        obj.multi_stim.data(ind).wf(roi_ind,:) = roi_value;
                    elseif metastats.stim{trial} == "N"
                        ind = length(obj.no_stim.data);
                        obj.no_stim.data(ind).wf(roi_ind,:) = roi_value;
                    end 

                    % Now write beh data
                    if metastats.beh{trial} == "H"
                        ind = length(obj.hit.data);
                        obj.hit.data(ind).wf(roi_ind,:) = roi_value;
                    elseif metastats.beh{trial} == "M"
                        ind = length(obj.miss.data);
                        obj.miss.data(ind).wf(roi_ind,:) = roi_value;
                    elseif metastats.beh{trial} == "FA"
                        ind = length(obj.false_alarm.data);
                        obj.false_alarm.data(ind).wf(roi_ind,:) = roi_value;
                    elseif metastats.beh{trial} == "CR"
                        ind = length(obj.correct_rejection.data);
                        obj.correct_rejection.data(ind).wf(roi_ind,:) = roi_value;
                    elseif metastats.beh{trial} == "EL"
                        ind = length(obj.early_lick.data);
                        obj.early_lick.data(ind).wf(roi_ind,:) = roi_value;
                    end   
                    

%                     if strcmpi(stim, 'V')
%                         obj.cti(4) = obj.cti(4)+1;
%                         obj.visual_stim.data(roi_ind,:,obj.cti(4)) = roi_value;
%                     elseif strcmpi(stim, 'S')
%                         obj.cti(5) = obj.cti(5)+1;
%                         obj.sensory_stim.data(roi_ind,:,obj.cti(5)) = roi_value;                  
%                     elseif strcmpi(stim, 'V+S')
%                         obj.cti(6) = obj.cti(6)+1;
%                         obj.multi_stim.data(roi_ind,:,obj.cti(6)) = roi_value;
%                     elseif strcmpi(stim, 'N')
%                         obj.cti(7) = obj.cti(7)+1;
%                         obj.no_stim.data(roi_ind,:,obj.cti(7)) = roi_value;
%                     else
%                         warning(strcat("Unknown stimulus type: ", stim));
%                     end
% 
%                     if strcmpi(beh, 'H')
%                         obj.cti(8) = obj.cti(8)+1;
%                         obj.hit.data(roi_ind,:,obj.cti(8)) = roi_value;
%                     elseif strcmpi(beh, 'M')
%                         obj.cti(9) = obj.cti(9)+1;
%                         obj.miss.data(roi_ind,:,obj.cti(9)) = roi_value;                  
%                     elseif strcmpi(beh, 'FA')
%                         obj.cti(10) = obj.cti(10)+1;
%                         obj.false_alarm.data(roi_ind,:,obj.cti(10)) = roi_value;
%                     elseif strcmpi(beh, 'CR')
%                         obj.cti(11) = obj.cti(11)+1;
%                         obj.correct_rejection.data(roi_ind,:,obj.cti(11)) = roi_value;
%                     elseif strcmpi(beh, 'EL')
%                         obj.cti(12) = obj.cti(12)+1;
%                         obj.early_lick.data(roi_ind,:,obj.cti(12)) = roi_value;
%                     else
%                         warning(strcat("Unknown response type: ", beh));
%                     end

                end
            end
        end
        
%         % Merge session data and index with the right class variable
%         if strcmpi(session_type, "visual_task")
%             obj.visual_task.data = session_data; obj.cti(1) = session_ind;
%             obj.visual_task.beh = session_beh;
%             obj.visual_task.stim = session_stim;
%             obj.visual_task.session = session_session;
%         elseif strcmpi(session_type, "sensory_task")
%             obj.sensory_task.data = session_data; obj.cti(2) = session_ind;
%             obj.sensory_task.beh = session_beh;
%             obj.sensory_task.stim = session_stim;
%             obj.sensory_task.session = session_session;
%         elseif strcmpi(session_type, "naive_task")
%             obj.naive_task.data = session_data; obj.cti(3) = session_ind;
%             obj.naive_task.beh = session_beh;
%             obj.naive_task.stim = session_stim;
%             obj.naive_task.session = session_session;
%         end
    end
end
    
end
