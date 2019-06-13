%% This function is deprecated and not needed anymore
function data = data_loader_one_condition(path, condition)
% Receives condition and path and loads the data

load(path);

switch condition
    case "hit"
        data = reshape([hit.data.wf], [15, 200, length(hit.data)]);
    case "miss"
        data = reshape([miss.data.wf], [15, 200, length(miss.data)]);
    case "false_alarm"
        data = reshape([false_alarm.data.wf], [15, 200, length(false_alarm.data)]);
    case "correct_rejection"
        data = reshape([correct_rejection.data.wf], [15, 200, length(correct_rejection.data)]);
    case "visual_stim"
        data = reshape([visual_stim.data.wf], [15, 200, length(visual_stim.data)]);
    case "sensory_stim"
        data = reshape([sensory_stim.data.wf], [15, 200, length(sensory_stim.data)]);
    case "no_stim"
        data = reshape([no_stim.data.wf], [15, 200, length(no_stim.data)]);
    case "multi_stim"
        data = reshape([multi_stim.data.wf], [15, 200, length(multi_stim.data)]);
    case "early_lick"
        data = reshape([early_lick.data.wf], [15, 200, length(early_lick.data)]);
    case "visual_task"
        data = reshape([visual_task.data.wf], [15, 200, length(visual_task.data)]);
    case "sensory_task"
        data = reshape([sensory_task.data.wf], [15, 200, length(sensory_task.data)]);
    case "naive_task"
        data = reshape([naive_task.data.wf], [15, 200, length(naive_task.data)]);
end 
        
end

