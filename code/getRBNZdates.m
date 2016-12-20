function [yr,qtr]=getRBNZdates(varargin)   
% PURPOSE: Extracts year and quarter in correct format to be used in name
% for storage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin==1    
    dat=varargin{1};
    currentQtr=['Q' int2str(rem(dat,1)*100) '-' int2str(floor(dat))];
else
    % get current quarter in format 27, eg	'QQ-YYYY' Q1-2001
    currentQtr=datestr(now,27);
end;

% quarter          
if strcmpi(currentQtr(:,1:2),'Q1')
    qtr='mar';
elseif strcmpi(currentQtr(:,1:2),'Q2')
    qtr='jun';
elseif strcmpi(currentQtr(:,1:2),'Q3')
    qtr='sep';
elseif strcmpi(currentQtr(:,1:2),'Q4')
    qtr='dec';
end;
% year
yr=currentQtr(:,end-1:end);