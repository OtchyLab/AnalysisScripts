function [records, ISI, list] = FindISI(cell)
%Function takes in the cell file and returns the record indices of those
%that have ISI's less than some threshhold (defined below)

thresh = 0.001; %Sets the max limit for which ISIs will be scanned

numrecs = length(cell);
records = zeros(numrecs,1);
ISI = -1*ones(numrecs,10);

for i=1:numrecs
    sortISI = sort(diff(cell(1,i).spikes));
    minISI= min(diff(cell(1,i).spikes));
    if minISI<=thresh
        records(i) = 1;
        ISI(i,:) = sortISI(1:10);
    end
end

list = find(records==1);

end