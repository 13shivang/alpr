function [r takethisbox]=connn(NR)

[Q,W]=hist(NR(:,4),8);
% Histogram of the y-dimension widths of all boxes.
%distributes NR in 20 bins
%Q is a row vector of frequency 
%and W is the row vector of all the mid
% points of bins.

q1=max(Q);
ind=find(Q==q1);
% Find indices from Q corresponding to frequency q1.
% q1 can be adjusted to th enumber of charcters in the numbr plate

if length(ind)==1  
    % If q1 boxes of interest are succesfully found record
    MP=W(ind);
    %  the midpoint of corresponding bin.
    binsize=W(2)-W(1);
    % Calculate the container size.
    container=[MP-(binsize/2) MP+(binsize/2)];
     % Calculating the complete container size
    [r takethisbox]=takeboxes(NR,container,2);
    %CHK=2 considers y-dimension width grouping
else
    r=[];
    takethisbox=[];
end
end