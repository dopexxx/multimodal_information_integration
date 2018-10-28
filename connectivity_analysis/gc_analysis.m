%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CONTEXT AND SOURCES:
% Class for granger causality analysis of Ca2+ data of multiple ROI of mice
%   behaving under different task conditions. Data similar to described in:
%    [1] "Context-dependent cortical integration of visual and somatosensory 
%           stimuli in behaving mice" M.Buchholz, Y.Sych, F.Helmchen, A.Ayaz
%  Utilized toolbox: Multivariate Granger Causality Analysis Toolbox:
% [2] L. Barnett, A. K. Seth, The MVGC Multivariate Granger Causality Toolbox: 
%   A New Approach to Granger-causal Inference, J. Neurosci. Methods (2014).
%
% FUNCTION:
% This class automatizes granger causality analysis of all sessions of a 
% given mouse that fulfill a given criterion (e.g. visual task). It expects
% to receive 2 criteria, so it can run the gc analysis for both sets of
% sessions separately and then perform statistical significance testing on 
% whether the connectivity differs from one criterion to another (e.g.
% visual vs. sensory task)
%
%   Jannis Born, October 2018

classdef gc_analysis

properties (Access = public)
    % Declare class variables
    mouse
    group_crit
    max_time_lag
end

properties (Access = private)
    % Internal class variables
   allowed_mice = ["5627rr", "5212r", "1110r", "1111lr", "1113rr", "2905l", "2907ll"];
   allowed_group_crits = ["S", "V", "N"];
   allowed_time_lag = [0, 1000];
   file_paths = ["I:\data\registered", "F:\data\registered", "/Users/jannisborn/Desktop/HIFO/multimodal_information_integration"];
   delimiter = '/'; % / for Mac, \ for Windows
   meta_file_root = "H:\data\behavior\";
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
    % (3) max_time_lag - {<double>, <intX>}. Specifies the maximal time 
    %           lag (maximum model order) for the GC analysis. How much
    %           time is maximally allowed for the signal to travel from
    %           region A to B. Give in ms. Defaults to 300ms, corres-
    %           ponding to 15 frames (SR: 20 Hz). max_time_lag is only
    %           an upper bound, ideal model order selected via BIC.

    % DEPRECATED
    % (4) padding - {<string>} from {"VALID", "SAME"}. [2] expects a
    %           matrix #ROI x frames_per_trial x #trials. In [1]
    %           frames_per_trial = 200 for all data, but #trials varies
    %           highly from session to session

    function obj = gc_analysis(varargin)

        % Error handling
        if (nargin < 2 || ~isstring(varargin{1}) || ~isstring(varargin{2})...
                || length(varargin{2})~= 2)
            error(['Please ensure the first arg is a STRING for the '...
                'mouse name and second arg is a STRING of length 2 '...
                'for the 2 groups to compare']);
        elseif nargin==3 && ~isnumeric(varargin(3))
            error("Please fed a numeric value for max_time_lag");
        elseif nargin > 3
            warning("Fourth and all later args are discarded");
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

        if nargin >= 3 && (varargin{3} < obj.allowed_time_lag(1) || ...
                varargin{3} > obj.allowed_time_lag(2))
            warning("Choose a maxmimum time lag [0, 1000]. Default is used instead")
            obj.max_time_lag = 300;
        elseif nargin == 2
            obj.max_time_lag = 300;
        else
            obj.max_time_lag = varargin{3};
        end

    end


    function assemble_data(obj)
    % This function assembles the data matrices. It searches for all
    % sessions of a given mouse and extracts the trials according to
    % the criteria defined in group_crit.

    for gc = obj.group_crit

        % For each hard-disk, go through all subfolders
        for file_path = obj.file_paths

            date_folders = obj.list_subfolders(file_path);

            % Open first level of folder (folder names are dates)
            for date_folder = date_folders

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
                            % TODO: If other grou criteria are allowed,
                            % this needs modification (optional execution
                            % only)
                            meta_path = strcat(obj.meta_file_root, date_id, ...
                                obj.delimiter, obj.mouse, obj.delimiter, ...
                                session_id);
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
                                session_data = obj.parse_data(metastats, ...
                                    session_folder);
                                % Append session data or do sth similar
                                
                            end
                                
                            

                        end
                    end
                end


            end


        end
    end


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
    
    function session_type = get_session_type(obj, metastats)
        % Receives the stats of a session and returns whether the session
        % was visual ("V"), somatosensory ("S") or naive ("N") --> ?.
        
        if sum(strcmpi(meta.stim,"S") & strcmpi(meta.beh,"H")) > 0 
            session_type = "V";
        elseif sum(strcmpi(meta.stim,"V") & strcmpi(meta.beh,"H")) > 0 
            session_type = "S";
        else
            error(strcat("Session type unclear for metastats ", metastats));
        end
        
    end
    
    function data = parse_data(obj, metastats, path)
        % Receives the path to a given session and the metastats of the
        % session. Loops over all trials, reads the CA2+ matrix, warps the 
        % data to the standard atlas, applies the ROI mask for each
        % preserved ROI, averages and returns a matrix of shape:
        % #ROI x 200 x num_trials
        
        42;
        
        
    end


end
    
end
