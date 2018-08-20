for i = 1:length(elements)
    numStarts = length(elements{i}.segFileStartTimes);
    numSegs = length(elements{i}.segType);
    if numStarts ~= numSegs
        display([num2str(i) '  --  ' keys{i}])
    end
end
