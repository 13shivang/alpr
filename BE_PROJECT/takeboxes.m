function [r takethisbox]=takeboxes(NR,container,chk)

takethisbox=[];
% Initialize the variable to an empty matrix.
for i=1:size(NR,1)
    if NR(i,(2*chk))>=container(1) && NR(i,(2*chk))<=container(2)
         %2*chk =4 gives y dimension width
         % If Bounding box is among the container plus tolerenc
        takethisbox=cat(1,takethisbox,NR(i,:));
        % Take that box and concatenate along first dimension
        % adding NR(properties)to take this box
    end
end
r=[];
for k=1:size(takethisbox,1)
    %till there are boxes in takethisbox (first dimensoin)
    var=find(takethisbox(k,1)==reshape(NR(:,1),1,[]));
    % Finding the indices of the interested boxes among NR
    if length(var)==1
        % since x-coordinate of the boxes will be unique.
        r=[r var];
    else% In case if x-coordinate is not unique
        for v=1:length(var)
            % then check which box fall under container condition. 
            M(v)=NR(var(v),(2*chk))>=container(1) && NR(var(v),(2*chk))<=container(2);
        end
        r=[r var];
    end
end
end