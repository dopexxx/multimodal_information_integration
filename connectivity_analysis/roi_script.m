%% Script to run the ROI extraction (Jannis Born, June 2019).
%
% This script is used to run the ROI extraction (see roi_extraction.m).

%% Specify meta variables
mouse = "5627rr";
disk = "Moritz_wide";

%%
obj = roi_extraction(mouse, disk);
obj.assemble_data();