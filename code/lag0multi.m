function out=lag0multi(x,p)


% INPUTS:
%
% X = vector or matrix of data
% P = Maximum lag to take


[R,C]=size(x);

%Take the first R-p rows of matrix x


out = zeros(R,p*C);

count = 1;
for ii = 1:p
 
    
    out(ii+1:end,count:count-1+C) = x(1:(R-ii),:);
    
    count = count+C;
    
end


%Revision History:
%v 1.01
%Replaced x1=x(1:(length(x)-p),:) with x1=x(1:(R-p),:); 
%Now should work correctly for matrices;