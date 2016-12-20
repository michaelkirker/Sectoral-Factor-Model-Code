function fu=convertFreqS(fi)
% PURPOSE: Convert from string to numeric

if ischar(fi)
    if strcmpi(fi,'a')
        fu=1;
    elseif strcmpi(fi,'q')
        fu=4;
    elseif strcmpi(fi,'m')
        fu=12;
    elseif strcmpi(fi,'d')
        fu=365;
    elseif strcmpi(fi,'b')
        fu=262;                
    end;
else
    fu=fi;
end;