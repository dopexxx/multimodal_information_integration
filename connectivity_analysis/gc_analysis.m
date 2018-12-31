function results = gc_analysis(data, sr, momax, save_path, brain_areas)
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
% This method automatizes granger causality analysis of all sessions of a 
% given mouse that fulfill a given criterion (e.g. visual task). It expects
% to receive
%`  (1) data - a 3D matrix of shape num_rois x num_obs x num_trials
%   (2) sr - a scalar indicating  the sampling rate (e.g. 20 for 20 Hz)
%   (3) momax - a scalar indicating the maximal time lag (e.g. 6 indicates
%           6 measurements, i.e. 300ms if sr = 20).
%   (4) save_path - a folder path were the data is saved
%   (5) brain_areas - a list of strings of length num_rois.
% It runs gc analysis and returns a struct with all necessary data for
% plotting.
%
%   Jannis Born, November 2018


%% Meta variables
regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
acmaxlags = 5000;   % maximum autocovariance lags (empty for automatic calculation)
tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.001;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

%% Start to work on the data

nrois = size(data,1);  % number of observations per trial
nobs = size(data,2);  % number of observations per trial
ntrials = size(data,3);     % number of trials

folder_name = strsplit(save_path,'/'); folder_name = folder_name(end-1);
mouse = strsplit(folder_name,'_'); mouse = mouse(1);
condition = strsplit(folder_name,'_'); condition = strjoin(condition(2:end-1));


results = struct();
results.data = data;
results.momax = momax;
results.sampling_rate = sr;

try
    %% Calculate information criteria
    % Calculate information criteria up to specified maximum model order.
    [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(data,momax,icregmode,false);

    % Plot information criteria.

    figure('units','normalized','outerposition',[0 0 1 1]); clf;
    plot_tsdata([AIC BIC]',{'AIC','BIC'},1/sr);
    title('Model order estimation');
    set(gcf,'color','w');
    savefig(strcat(save_path,"AIC_BIC_model_order_estimation.fig"));
    clf;

    fprintf('\nbest model order (AIC) = %d\n',moAIC);
    fprintf('best model order (BIC) = %d\n',moBIC);

    % Select model order.
    if strcmpi(morder,'AIC')
        morder = moAIC;
        fprintf('using AIC best model order = %d\n',morder);
    elseif strcmpi(morder,'BIC')
        morder = moBIC;
        fprintf('using BIC best model order = %d\n',morder);
    end
    results.morder = morder;

    %% VAR model estimation (<mvgc_schema.html#3 |A2|>)

    % Estimate VAR model of selected order from data.

    ptic('*** tsdata_to_var... ');
    [results.A,results.SIG] = tsdata_to_var(data,morder,regmode);
    ptoc;

    % Check for failed regression

    assert(~isbad(results.A),'VAR estimation failed');

    % NOTE: at this point we have a model and are finished with the data! - all
    % subsequent calculations work from the estimated VAR parameters A and SIG.


    %% Autocovariance calculation (<mvgc_schema.html#3 |A5|>)

    % The autocovariance sequence drives many Granger causality calculations (see
    % next section). Now we calculate the autocovariance sequence G according to the
    % VAR model, to as many lags as it takes to decay to below the numerical
    % tolerance level, or to acmaxlags lags if specified (i.e. non-empty).

    ptic('*** var_to_autocov... ');
    [results.G,results.info] = var_to_autocov(results.A,results.SIG,acmaxlags);
    ptoc;

    % The above routine does a LOT of error checking and issues useful diagnostics.
    % If there are problems with your data (e.g. non-stationarity, colinearity,
    % etc.) there's a good chance it'll show up at this point - and the diagnostics
    % may supply useful information as to what went wrong. It is thus essential to
    % report and check for errors here.

    var_info(results.info,true); % report results (and bail out on error)

    %% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

    % Calculate time-domain pairwise-conditional causalities - this just requires
    % the autocovariance sequence.

    % This takes time
    ptic('*** autocov_to_pwcgc... ');
    F = autocov_to_pwcgc(results.G);
    results.F = F;
    ptoc;

    % Check for failed GC calculation

    assert(~isbad(results.F,false),'GC calculation failed');

    %% Significance test using theoretical null distribution, adjusting for multiple
    % hypotheses.
    pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nrois-2,tstat);
    sig  = significance(pval,alpha,mhtc);

    % Plot time-domain results
    % (1) causal graph
    figure('units','normalized','outerposition',[0 0 1 1]); clf;
    set(gcf,'color','w');
    plot_pw(F, [], brain_areas);
    title(strcat('Pairwise-conditional GC, mouse = ', mouse, ', condition = ', condition));
    savefig(strcat(save_path,"Pairwise_gc"));
    % (2) p-values  
    figure('units','normalized','outerposition',[0 0 1 1]); clf;
    set(gcf,'color','w');
    plot_pw(pval, [], brain_areas);
    title(strcat('p-values, mouse = ', mouse, ', condition = ', condition));
    savefig(strcat(save_path,"p-values"));

    % (3) significance.
    figure('units','normalized','outerposition',[0 0 1 1]); clf;
    set(gcf,'color','w');
    plot_pw(sig, [], brain_areas);
    title(strcat('Significant at p = ', num2str(alpha),', mouse = ', mouse, ', condition = ', condition));
    savefig(strcat(save_path,"significance"));

    % Joint plot
    figure('units','normalized','outerposition',[0 0 1 1]); clf; set(gcf,'color','w');
    subplot(1,3,1);
    plot_pw(F, [], brain_areas);
    title(strcat('Pairwise-conditional GC, mouse = ', mouse, ', condition = ', condition));
    subplot(1,3,2);
    plot_pw(pval, [], brain_areas);
    title(strcat('p-values, mouse = ', mouse, ', condition = ', condition));
    subplot(1,3,3);
    plot_pw(sig,[], brain_areas);
    title(strcat('Significant at p = ', num2str(alpha),', mouse = ', mouse, ', condition = ', condition));
    savefig(strcat(save_path,"connectivity"));

    % For good measure we calculate Seth's causal density (cd) measure - the mean
    % pairwise-conditional causality. We don't have a theoretical sampling
    % distribution for this.

    cd = mean(F(~isnan(F)));

    fprintf('\ncausal density = %f\n',cd);
    results.causal_density = cd;

    % Save the local variables
    save(strcat(save_path,'data'),'results');
    close all;
catch
    warning("Error in gc_analysis.m")
end
    
end

