%% CLEAN DATA
%
% This code reads in the raw data, cleans it, transforms it, and saves it
% down into a mat file for the estimation function to use


function cleandata(info)
% Reads in raw data from excel spreadsheet. Cleans and transforms data and
% outputs it into a structure ready to use in the estimation.


% estimation start details
begyr   = floor(info.startEstimation);
begqtr = round(rem(info.startEstimation,1)*100);

% estimation end details
endyr   = floor(info.endEstimation);
endqtr = round(rem(info.endEstimation,1)*100);

timeVector = latttt(begyr,endyr,begqtr,endqtr,'freq',4);








%% ===================================================================== %%
%   READ IN DATA
%  =====================================================================  %

% Read in the price indices

% Load CPI index numbers (stored in Sheet1)
[rawdata, text]=xlsread([info.inputDir '\price_indices.xls'],'Sheet1');

names = text(1,2:end); % Names of the series in the spreadsheet


% Load in the series loadinga for each factor
[fac_loadings]=xlsread([info.inputDir '\factor_loading_classification.xls'],'loadings'); %Tradables and non-tradable dummys



%% ===================================================================== %%
%                           2) Transform data
%  =====================================================================  %


% Raw data may be longer than the estimation period. So trim initial
% observations based on length of timeVector (the sample to be used in
% estimation)
rawdata=rawdata(end-length(timeVector)+1:end,:);

% Transform raw data indices to APC and store for potential analysis later
data.rawdata = (rawdata(5:end,:)-rawdata(1:end-4,:))./rawdata(1:end-4,:)*100;

%-------------------------------------------------------------------------%
%   Adjust data: Interpolation
%-------------------------------------------------------------------------%
if info.Interp % If user requested interpolation of particular date(s):
    
    
    [~,~,Adates] = intersect(info.interpdates,timeVector); % See if adjustment date(s) are in the sample
    
    if ~isempty(Adates)
        for ii = 1:length(Adates) % make adjustment for each date
            
            qpc1 = rawdata(Adates(ii)-1,:)./rawdata(Adates(ii)-2,:)-1; % QPC of quarter before the interpolation
            qpc2 = rawdata(Adates(ii)+1,:)./rawdata(Adates(ii),:)-1; % QPC of the quarter after the interpolation
            
            newindx  = rawdata(Adates(ii)-1,:).*(1+(qpc1+qpc2)/2); % Interpolated indice value
            oldindx = rawdata(Adates(ii),:);
            
            % Since we have interpolated a new level for a particular
            % quarter, we need to shift all following quarters by the same
            % amount.
            rawdata(Adates(ii):end,:) = rawdata(Adates(ii):end,:)-repmat((oldindx-newindx),length(timeVector)-Adates(ii)+1,1);
        end
    end
    
    % Store a list of the dates that were interpolated
    data.interpdates = timeVector(Adates);
end

%-------------------------------------------------------------------------%
% Convert data to APC
%-------------------------------------------------------------------------%
infl = (rawdata(5:end,:)-rawdata(1:end-4,:))./rawdata(1:end-4,:)*100;

ann_timeVector = timeVector(5:end); % trim time vector to match APC data length


%-------------------------------------------------------------------------%
% Adjust data: remove GST
%-------------------------------------------------------------------------%
if info.GSTadjust
    
    [~,~,Adates] = intersect(info.adjustdates,timeVector(5:end)); % See if GST period is in the sample
    
    if ~isempty(Adates)
        for ii = 1:length(Adates)
            
            x_max = min(Adates+3,length(timeVector(5:end))); % because the model runs off annual inflation, we need to adjust the next 3 quarters worth of data as well
            dummymoments  = xlsread([info.inputDir '/GSTadjustments.xls'],num2str(ann_timeVector(Adates)));
            
            infl(Adates(ii):x_max,:)= infl(Adates(ii):x_max,:) - repmat(dummymoments',length(Adates(ii):x_max),1);
        end
    end
    
    % store dates GST adjustment was made.
    data.GSTdates = ann_timeVector(Adates);
end



%-------------------------------------------------------------------------%
%   Remove series not spanning entire sample length
%-------------------------------------------------------------------------%

toremove = find(isnan(sum(infl))==1); % Fins which series have missing data

% Drop series
infl(:,toremove)    = [];
names(toremove)     = [];
fac_loadings(toremove,:)  = [];
data.rawdata(:,toremove) = [];


names = lower(names);

% Find the series the relate to CPI (pcpcis), Nontradables (pnt), and
% Tradables (ptr)

pcpisindx = strmatch('pcpis',names,'exact');
pntindx = strmatch('pnt',names,'exact');
ptrindx = strmatch('ptr',names,'exact');

data.keyindx = [pcpisindx ptrindx pntindx]; % location of key series used to identify model


%-------------------------------------------------------------------------%
%                   Standardise the data
%-------------------------------------------------------------------------%

data.sermean = mean(infl);
data.serstd  = std(infl);

data.infl = standardise(infl);

data.names = names;
data.fac_loadings = fac_loadings;
data.dates = ann_timeVector;

%=========================================================================%
%   Save cleaned data structure
%=========================================================================%

[yr,qtr]=getRBNZdates(info.endEstimation);

save([info.tempDir '\' 'data_' yr qtr '.mat'],'data');



