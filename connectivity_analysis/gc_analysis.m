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

properties (Access = public)
    % Declare class variables
    root = "H:\Jannis\granger_causality_data\";
    delim = "\";
end



 % Error handling
        if (nargin < 2 || ~isstring(varargin{1}) || ~isstring(varargin{2})...
                || length(varargin{2})~= 2)
            error(['Please ensure the first arg is a STRING for the '...
                'mouse name and second arg is a STRING of length 2 '...
                'for the 2 groups to compare']);
            
            
if         elseif nargin==3 && ~isnumeric(varargin(3))
            error("Please fed a numeric value for max_time_lag");
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

        end
end