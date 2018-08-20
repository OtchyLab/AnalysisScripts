function diffVec = timeDifferences(filename)
    % function to get a vector of time differences for the data readings in a
    % .mov file
    % takes the name of a csv file that corresponds to a .mov file
    % read in csv file
    data = csv2cell(filename, 'fromfile');
    % first 13 rows do not contain numerical data; go 14-end of col 1
    [r, ~] = size(data); 
    timeStamps = data(14:r,1);
    diffVec = zeros(1, length(timeStamps)-1);
    for i = 1:length(timeStamps)-1
        diffVec(i) = str2num(timeStamps{i+1}) - str2num(timeStamps{i});
    end
end