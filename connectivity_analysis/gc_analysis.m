%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CONTEXT AND SOURCES:
% Class for granger causality analysis of Ca2+ data of multiple ROI of mice
%   behaving under different task conditions. 
%  Utilized toolbox: Multivariate Granger Causality Analysis Toolbox:
%   [1] L. Barnett, A. K. Seth, The MVGC Multivariate Granger Causality Toolbox: 
%       A New Approach to Granger-causal Inference, J. Neurosci. Methods (2014).
%
% FUNCTION:
% This class automatizes granger causality analysis of all sessions of a 
% given mouse that fulfill a given criterion (e.g. visual task). It expects
% to receive the paths to 2 '.m' files that contain data for two criteria
% (e.g. visual task vs sensory task). It runs gc analysis for both sets of
% observations separately and then perform statistical significance testing 
% whether the connectivity differs from one criterion to another.
%
%   Jannis Born, November 2018

classdef gc_analysis

properties (Access = public)
    % Declare class variables
    path_1
    path_2
    max_time_lag
end

    
end
