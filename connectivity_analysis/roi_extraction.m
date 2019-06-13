%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CONTEXT AND SOURCES:
% Class for data gathering of Ca2+ widefield data of multiple ROI 
%   of mice behaving under different task conditions. 
%   Data similar to:
%   [1] "Context-dependent cortical integration of visual and somatosensory 
%           stimuli in behaving mice" M.Buchholz, Y.Sych, F.Helmchen, A.Ayaz
% 
% FUNCTION:
% This is a preparatory of the gc_analysis class. Extracted ROIs can be
% used for all types of analysis (correlation etc.)
% INPUTS:
%   (1) mouse -  {<string>}, specifying name of the mouse. Choose from
%       {"5627rr", "5212r", "1110r", "1111lr", "1113rr", "2905l", "2907ll"}.
%       Mouse ID for which analysis is performed.
%   (2) disk - {<string>, <char>}, specifying the disk to search on. 
%
% SAVES:
%   For each of the following groups (see below), one struct is saved in
%       the specified path.
%
%   Each file has shape #ROI x samples_per_trial x num_trial 
%
%   NOTE:
%      Sessions where the mouse was unreactive (always or never licked) 
%      are discarded.
%   
%   NOTE 2:
%       Please CHECK whether your input data was already co-registered. If
%       YES, you need to comment out the 'warper' variable, if NOT, you the
%       object for the affine transformation has to be in the paths given
%       to the variable 'affine'.
%
%   Jannis Born, October 2018

classdef roi_extraction

properties (Access = public)
    % Declare class variables
    mouse
    disk
    sensory_task; visual_task; naive_task; 
end

properties (Access = private)
   
    % List of allowed_mice and disks is to be extended.
   allowed_mice = ["5627rr", "5212r", "1110r", "1111lr", "1113rr", "2905l", "2907ll", "2906r"];
   allowed_disks = ["Moritz_wide", "Moritz_img2", "F", "I"];
   

   % Path variables
   delimiter = '/'; % / for Mac, \ for Windows
   data_root_pre = "/Volumes/"
   data_root_post = "/data/registered";
   
   % Points to the disk/path with the behavioral result files
   meta_file_root = "/Volumes/Moritz_beh/data/behavior/";
   save_path = "/Volumes/Moritz_beh/Jannis/widefield_for_network/";
  
   % Load the ROIs and the affine transformation matrices for the warping
   % (one per mouse, needs to be done prior to analysis! Morio did that.)
   affine = load('/Users/jannisborn/Desktop/HIFO/clean/resources/transformation_matrices');
   rois = load('/Users/jannisborn/Desktop/HIFO/clean/resources/rois');
   
   % Check that this list matches the ROIs loaded above
   brain_areas = ["PL","ACC","M2","M1","S1b","S1mn","S1bf","S2","Au", ...
       "ASA","V2L","RL","V2A","V1","RS"];

   % Widefield image data resolition
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

    function obj = roi_extraction(varargin)

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
        
        if any(contains(obj.allowed_disks, varargin{2}))
            obj.disk = varargin{2};
        else
             error(strcat("Unknown disk name given (", varargin{2}, ...
                ") see help for details.")); 
        end

        %affine2d class to warp to standard atlas
        obj.warper = eval(strcat('obj.affine.transform.mouse',obj.mouse)); 
        
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
      
        obj.sensory_task = struct();
        obj.sensory_task.mouse = obj.mouse;
        obj.sensory_task.rois = obj.brain_areas;
        obj.sensory_task.data = struct();
        
        obj.visual_task = struct();
        obj.visual_task.mouse = obj.mouse;
        obj.visual_task.rois = obj.brain_areas;
        obj.visual_task.data = struct();
       
        obj.naive_task = struct();
        obj.naive_task.mouse = obj.mouse;
        obj.naive_task.rois = obj.brain_areas;
        obj.naive_task.data = struct();
          
    end
    
    function save_all(obj)
    
    % save all variables
    visual_task = obj.visual_task;
    sensory_task = obj.sensory_task;
    naive_task = obj.naive_task;
    save(strcat(obj.save_path,'gca_visual_task_', obj.mouse, '_disk_', obj.disk),'visual_task','-v7.3');
    save(strcat(obj.save_path,'gca_sensory_task_', obj.mouse, '_disk_', obj.disk),'sensory_task','-v7.3');
    save(strcat(obj.save_path,'gca_naive_task_', obj.mouse, '_disk_', obj.disk),'naive_task','-v7.3');
    
    end
    
    function obj = parse_data(obj, path, session_type, metastats)
        % Receives the path to a given session and the metastats of the
        % session. Loops over all trials, reads the CA2+ matrix, warps the 
        % data to the standard atlas, applies the ROI mask for each
        % preserved ROI, averages and returns a matrix 'data' of shape:
        % #ROI x samples_per_trial x num_trials
        registration = load(strcat(path,obj.delimiter,'registration'));
   
 
        % Error handling
        if ~(registration.info.trials_obj==length(metastats.beh))
            warning(strcat("Found CA data for ", num2str(registration.info.trials_obj), ...
                " trials, but behavioral data for ", num2str(length(metastats.beh)), ...
                ". Will neglect overhead trials."));
        end
        
        tic
        % Load imaging data trial per trial
        for trial = 1:min(registration.info.trials_obj,length(metastats.beh))
            error = 0;
            if mod(trial, 1000) == 0
                toc
                disp(['Currently processing trial ', num2str(trial),...
                    '/', num2str(length(metastats.beh))]);
                tic
            end
            try 
                % Load trial data (use matfile since it is faster than load)
                trial_data = matfile(strcat(path,obj.delimiter,'dFF_t',num2str(trial)));

                data = trial_data.dFF;
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
                
                % Comment this line out if warping/coregistration was done already
                data = imwarp(data, obj.warper,'OutputView', obj.cutter);
                data = reshape(data, [size(data,1)*...
                    size(data,2), size(data,3)]);

                % FIRST: Allocate space and write behavioral/stimulus data. 

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

                % Loop over ROIs and save values.
                % Save average response of all rois for all frames of all trials.
                for roi_ind = 1:length(obj.brain_areas)
                    
                    % The below line is needed if ROIs refers to the old
                    % set of 33 ROIs
                    %area_mask = eval(strcat('obj.rois.ROIs.',obj.brain_areas(roi_ind),'.maskCircle'));
                    
                    % This code is for the new 15 ROIs
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

                    
                end
            end
        end
    end
end
    
end
