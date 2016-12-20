function core = createcore(results,data,serindx)
%% ===================================================================== %%
%%           CREATE CORE INFLATION ESTIMATE FROM DFM RESULTS             %%
%% ===================================================================== %%
%
%   Created by: Michael Kirker - Modelling Team, RBNZ
%               michael.kirker@rbnz.govt.nz
%   Date last edited: 18/06/2011
%
% This functions takes the estimated DFM results structure, and constructs
% the core inflation estimate by scaling the factors by the appropriate
% factor loadings and series mean and standard deviation.
%
% INPUTS:
% 
% [RESULTS]: Structure containing the results from DFM estimation
%
% [DATA]: Structure containing the details of the data used to estimate the
% model. Must contain DATA.NAMES (an NNx1 cell of series names),
% DATA.SERMEAN (an NNx1 vector of mean inflation rates for each series),
% and DATA.SERSTD (an NNx1 vector of standard deviations for each series).
%
% [SERINDX]: (optional) - index number of series you want a core inflation
% measure for. If none is defined, the default is for the first series in
% the panel





% If no specific series are defined in SERINDX, compute a core measure of all by default
if nargin < 2
    serindx = 1;
    disp('You did not specify which core series. Therefore I will compute core for the first series.')
end




core = struct();



F = results.floads(serindx,:,:);

N = size(results.factors,1);

for ii = 1:length(serindx)
   
    eval(['core.' data.names{serindx(ii)} ' = data.sermean(serindx(ii)) + data.serstd(serindx(ii)).*sum(results.factors.*repmat(F(ii,:,:),N,1),3);']);
    
end