function options=validateInput(default,varargin)
% PURPOSE: Validate input. General function. Can be used by all functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(default)
    error('validateInput:err','default input needs to be passed as a cell');
end;

% extract defualt input 
defaultNames=default(1:3:end);
defaultValues=default(2:3:end);
validate=default(3:3:end);

default=cell2struct(defaultValues,defaultNames,2);
options=default;

% overwrite input with varargin values if present 
if nargin>1    
    optInput=varargin{1};            
    % NOTE: varargin can be a structure. Sometimes, if varargin is passed
    % from one function to a other function it can end up as a structure
    % inside a cell. The try catch block bellow fix this. If the structure
    % is inside a cell inside a cell, the function can not help you!    
    if isstruct(optInput) 
        optInNames=fieldnames(optInput);
        optInValues=struct2cell(optInput);                
    else    
        optInNames=optInput(1:2:end);        
        optInValues=optInput(2:2:end);
    end;    
    if ~iscellstr(optInNames)        
        try
            optInNames=fieldnames(optInput{:});
            optInValues=struct2cell(optInput{:});                                                    
        catch ME
            disp('validateInput:err','varargin needs to be a cell array with every second argument as a string')            
            rethrow(ME)
        end;
    end;            
    
    [matchingNames,ia,ib]=intersect(defaultNames,optInNames);
    
    % to be used for validation
    notValid=false(1,length(ia));validateErrorMsg='Following error(s) found in input: ';
    for mn=1:length(ia)        
        %validate input
        x=validate{ia(mn)};
        notValid(mn)=~x(optInValues{ib(mn)});       
        if ~notValid(mn)
            options.(defaultNames{ia(mn)})=optInValues{ib(mn)};    
        else
            validateErrorMsg=sprintf('%s\n%s',validateErrorMsg,['Input: ' defaultNames{ia(mn)} ', should have passed: ' func2str(x)]);
        end;
    end;    
    
    % print error message
    if any(notValid)  
        error('validateInput:err',validateErrorMsg);
    end;
end;



