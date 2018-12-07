%% GC ANALYSIS SCRIPT
%% Define task options
mouse = "5212r";

condition1 = "visual_task";
condition2 = "sensory_task";

%%
root = "H:\Jannis\granger_causality_data\";
delim = "\";
disk = "I";

brain_areas = ["V1 rostral", "V1 middle lateral", "V1 middle medial", ...
       "V1 caudal lateral", "V1 caudal medial", "V rostral", "V anterolateral", ...
       "V lateral", "V posteromedial", "V anterior", "Retrosplenial caudal", ...
       "Retrosplenial rostral", "Retrosplenial inter", "Aud", "S1 rostral", ...
       "S1 lateral", "S1 medial", "S1 caudal", "S trunk", "S hindlimb", ...
       "S forelimb", "S nose", "S mouth", "S supplementary", "M1 lateral", ...
       "M1 intermediate", "M1 medial", "M2 rostrolateral", "M2 rostromedial", ...
       "M2 intermiedate", "M2 caudal", "PFC prelimbic", "PFC ACC"];

p1 = strcat(root,"gca_",condition1,"_",mouse,"_disk_",disk);
p2 = strcat(root,"gca_",condition2,"_",mouse,"_disk_",disk);

load(p1); % visual_task.data
load(p2); % sensory_task.data

% Only investigate matrix 1
X = eval(strcat(condition2,".data"));

T = reshape(X, [size(X,1),size(X,2)*size(X,3)]);
missing_rois = find(~any(T,2));
disp(strcat('Exclude ROI numbers = ', num2str(missing_rois), ' due to missing data.'));

X(missing_rois,:,:) = [];
brain_areas(missing_rois) = [];

%% GCA FUNCTION
ntrials   = size(X,3);     % number of trials
nobs      = size(X,2);  % number of observations per trial
nrois      = size(X,1);  % number of observations per trial

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 6;     % maximum model order for model order estimation

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.001;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = 200;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

seed      = 0;      % random seed (0 for unseeded)

%%

% Calculate information criteria up to specified maximum model order.

ptic('\n*** tsdata_to_infocrit\n');
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
ptoc('*** tsdata_to_infocrit took ');

% Plot information criteria.

figure(1); clf;
plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
title('Model order estimation');

%amo = size(AT,3); % actual model order

fprintf('\nbest model order (AIC) = %d\n',moAIC);
fprintf('best model order (BIC) = %d\n',moBIC);
%fprintf('actual model order     = %d\n',amo);

% Select model order.

if     strcmpi(morder,'actual')
    morder = amo;
    fprintf('\nusing actual model order = %d\n',morder);
elseif strcmpi(morder,'AIC')
    morder = moAIC;
    fprintf('\nusing AIC best model order = %d\n',morder);
elseif strcmpi(morder,'BIC')
    morder = moBIC;
    fprintf('\nusing BIC best model order = %d\n',morder);
else
    fprintf('\nusing specified model order = %d\n',morder);
end

%% VAR model estimation (<mvgc_schema.html#3 |A2|>)

% Estimate VAR model of selected order from data.

ptic('\n*** tsdata_to_var... ');
[A,SIG] = tsdata_to_var(X,morder,regmode);
ptoc;

% Check for failed regression

assert(~isbad(A),'VAR estimation failed');

% NOTE: at this point we have a model and are finished with the data! - all
% subsequent calculations work from the estimated VAR parameters A and SIG.

%% Autocovariance calculation (<mvgc_schema.html#3 |A5|>)

% The autocovariance sequence drives many Granger causality calculations (see
% next section). Now we calculate the autocovariance sequence G according to the
% VAR model, to as many lags as it takes to decay to below the numerical
% tolerance level, or to acmaxlags lags if specified (i.e. non-empty).

ptic('*** var_to_autocov... ');
[G,info] = var_to_autocov(A,SIG,acmaxlags);
ptoc;

% The above routine does a LOT of error checking and issues useful diagnostics.
% If there are problems with your data (e.g. non-stationarity, colinearity,
% etc.) there's a good chance it'll show up at this point - and the diagnostics
% may supply useful information as to what went wrong. It is thus essential to
% report and check for errors here.

var_info(info,true); % report results (and bail out on error)

%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

% Calculate time-domain pairwise-conditional causalities - this just requires
% the autocovariance sequence.

% This takes time
ptic('*** autocov_to_pwcgc... ');
F = autocov_to_pwcgc(G);
ptoc;

% Check for failed GC calculation

assert(~isbad(F,false),'GC calculation failed');
%%
% Significance test using theoretical null distribution, adjusting for multiple
% hypotheses.
alpha = 0.001;
pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nrois-2,tstat); % take careful note of arguments!
sig  = significance(pval,alpha,mhtc);

% Plot time-domain causal graph, p-values and significance.

figure('units','normalized','outerposition',[0 0 1 1]); clf;
plot_pw(F, [], brain_areas);
title('Pairwise-conditional GC');

% For good measure we calculate Seth's causal density (cd) measure - the mean
% pairwise-conditional causality. We don't have a theoretical sampling
% distribution for this.

cd = mean(F(~isnan(F)));

fprintf('\ncausal density = %f\n',cd);
