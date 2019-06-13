function out = data_loader_task_specific(path, task, stim, beh)
% Receives condition and path and loads the data.
%
% ARGS:
%   PATH    - Path to the data.
%   TASK    - Choose from "visual_task", "sensory_task" and "naive_task"
%   STIM    - Defaults to "A" (all), options are "V", "S", "V+S", "N"
%   BEH     - Defaults to "A" (all), options are "H", "M", "CR", "EL", "FA"


if stim == "A" && beh == "A"
    warning("All trials will be extracted (independent of stimulus and behavior).")
end

load(path);

switch task
    case "visual_task"
        data = visual_task.data;
    case "sensory_task"
        data = sensory_task.data;
    case "naive_task"
        data = naive_task.data;
end 


keep_inds = [];
for trial = 1:length(data)
    if stim == "A" || data(trial).stim == stim
        if beh == "A" || data(trial).beh == beh
            keep_inds = [keep_inds, trial];
        end
    end
end

d = [data(keep_inds).wf];
out = reshape(d, [size(d,1), 200, length(keep_inds)]);
