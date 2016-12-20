%% SECTORAL FACTOR MODEL
%
% The code replicated the model in Kirker DP2010/13 "What drives core
% inflation? A dynamic factor model analysis of tradable and nontradable
% prices". This version contains fixes of typos present in the DP version
% of the model.
%
% This file runs all the code. Each section can be run indepedently
% provided the user has previously run the previous sections. Therefore,
% the user can comment out the estimation section, and re-run the graph
% code without having to re-estimate the entire model.


%% Housekeeping

home;
close all;
clear;


%% User input/settings

info = struct(); % Structure to contain user settings and other info.

%----------------------------------------------------------------------%
%                               DATES
%----------------------------------------------------------------------%
info.startEstimation = 1992.01;   % Start of data (YYYY.QQ)
info.endEstimation   = 2011.01;   % Final period of data (YYYY.QQ)
%----------------------------------------------------------------------%

info.progressbar = true; % Do you want a progress bar displayed?


%-------------------------------------------------------------------------%
%                   BAYESIAN ESTIMATION INPUTS
%-------------------------------------------------------------------------%
info.draws = 2000;  % number of draws
info.burn  = 1500;  % number of draws to burn

info.numfact = [1 1];   % The first element is the number of tradable factors
% to include in the model, the second element is the
% number of nontradable factors to include in the
% model.
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
%                           DATA ADJUSTMENT
%-------------------------------------------------------------------------%
% Changes in taxes (expecially broad-based taxes like GST) can dramatically
% affect the estimation of the model. These one-off price adjustments can
% bias the historical correlations between series the the factor model
% relies on.
%
% Below are two ways to correct for tax changes (if you wish). The first is
% to manually specify the size of the tax change (in QPC) that should be
% removed from each series as a result of the tax. The second is to
% interpolate the data for the quarter rather than use the actual.


% Take out GST effect
info.GSTadjust = true;
info.adjustdates = [2010.04];
% Create a list of dates that you want to interpolate. For each date, create
% a worksheet in the "GSTadjustments.xls" file, and give the QPC adjustment
% sizes to be removed for each series.


% Interpolate data for quarter
info.Interp = false;
info.interpdates = [2000.03];
% Interpolate a period where the size of GST is unkown by averaging the QPC
% of the quarters either side
%-------------------------------------------------------------------------%



%% Code Paths

restoredefaultpath;

info.currDir = cd;

addpath([info.currDir '/code']);

info.inputDir = [info.currDir '/input'];
info.tempDir = [info.currDir '/temp'];
info.outputDir = [info.currDir '/output'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
disp(' ');
disp(['Now running for: num of factors = ' num2str(info.numfact) ', date = ' num2str(info.endEstimation)]);
disp(' ');

[yr,qtr]=getRBNZdates(info.endEstimation);



%=========================================================================%
%     Clean Data
%=========================================================================%
% Read in raw data from spreadsheet, clean and transform the data, and
% store in mat file ready for the estimation code to use

cleandata(info)





%=========================================================================%
%       Estimate core inflation model
%=========================================================================%

% Read in cleaned data file
load([info.tempDir '/data_' yr qtr '.mat']);

% Estimate the model
[results] = DFMcore(data.infl,info.numfact,data.fac_loadings,info.draws,info.burn,data.keyindx,info.progressbar,1);


% Save down estimation results
results.dataused = ['data_' yr qtr '.mat'];

save([info.outputDir '/estimation_results/estimation_fac' num2str(info.numfact(1)) num2str(info.numfact(2)) '_' yr qtr '.mat'],'results')






%=========================================================================%
%       Construct core inflation measures from estimation
%=========================================================================%
% Load estimation results, and compute core inflation estimates for CPI, TR
% and NT inflation. Core inflation is computed for each saved draw of the
% Bayesian estimation
% {
load([info.tempDir '/data_'  yr qtr '.mat']);
load([info.outputDir '/estimation_results/estimation_fac' num2str(info.numfact(1)) num2str(info.numfact(2)) '_' yr qtr '.mat']);

T = length(data.dates);

core = createcore(results,data,data.keyindx);

core.dates = data.dates;
core.names = data.names(data.keyindx);
core.infl = repmat(data.sermean(data.keyindx),T,1) +  repmat(data.serstd(data.keyindx),T,1).*data.infl(:,data.keyindx);
core.rawdata = data.rawdata(:,data.keyindx);

save([info.outputDir '/core_inflation/core_fac' num2str(info.numfact(1)) num2str(info.numfact(2)) '_' yr qtr '.mat'],'core');
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below are a few simple outputs from the model.


%=========================================================================%
% Plot core inflation against headline inflation
%=========================================================================%

load([info.outputDir '/core_inflation/core_fac' num2str(info.numfact(1)) num2str(info.numfact(2)) '_' yr qtr '.mat']);


year_vec = floor(core.dates);
quarter_vec = (core.dates-year_vec)*100;

date_vec = year_vec + (quarter_vec-1)*0.25;

% CPI 
figure('name', 'CPI Core inflation vs Headline')
[p1 hh] = plotx2(date_vec(3:end)',prctile(core.pcpis ,[50 10 90],2));
hold on
pinf = plot(date_vec(3:end),core.rawdata(3:end,1),'b','LineWidth',1.3);


% Tradable
figure('name', 'Tradable Core inflation vs Headline Tradables')
[p1 hh] = plotx2(date_vec(3:end)',prctile(core.ptr ,[50 10 90],2));
hold on
pinf = plot(date_vec(3:end),core.rawdata(3:end,2),'b','LineWidth',1.3);


% Non-Tradable
figure('name', 'Nontradable Core inflation vs Headline Nontradables')
[p1 hh] = plotx2(date_vec(3:end)',prctile(core.pnt ,[50 10 90],2));
hold on
pinf = plot(date_vec(3:end),core.rawdata(3:end,3),'b','LineWidth',1.3);


%=====================================================================%
% Plot factor contribution
%=====================================================================%
% Plot the contribution of each factor type to core CPI inflation
load([info.outputDir '/estimation_results/estimation_fac' num2str(info.numfact(1)) num2str(info.numfact(2)) '_' yr qtr '.mat']);
plot_contrib_graphs(info,results,core); 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time2complete=toc;
time_hr  = floor(time2complete/3600);
time_min = floor((time2complete - time_hr*3600)/60);
time_sec = time2complete-time_hr*3600-time_min*60;

disp(['Code complete. Time taken to run = ' num2str(time_hr) ' hours, ' num2str(time_min) ' mins, ' num2str(time_sec) ' seconds.'])

