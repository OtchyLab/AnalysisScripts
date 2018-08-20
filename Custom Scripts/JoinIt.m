function JoinIt(vars, files, output)
%Function takes an ordered matrix of .mat files and concatenates the
%variables specified by the matrix vars.  If 'all' is specified, all 
%variable present in the first .mat file will be used. The concatenated
%variables are saved at location given by output.

if length(files)<2
    error('Must have more than one .mat file specified')
end

if length(vars)<2
    error('For simplicity, you must have more than one variable specified')
end

%Setup the list of variables to be processed
if size(vars,1) == 1 && strcmp(vars, 'all')
    active_vars = whos('-file', files{1});
    for i = 1:size(active_vars,1)
        buffer{i} = active_vars(i).name;
    end
    active_vars = buffer;
else
    active_vars = vars;
end

%active_vars = {active_vars;'end'};

%Determine number of variables to be processed
num_vars = size(active_vars,2);

%Suck out the specified vars from the first file into the *_new buffer
% if num_vars>1
     load(files{1}, active_vars{:})
% else
%     load(files{1}, active_vars)
% end

for i = 1:num_vars
    eval([active_vars{i} '_new=' active_vars{i} ';']);
end

%Iteratively load .mat files, append to the end of *_new
for i = 2:length(files)
    load(files{i}, active_vars{:});
    for j = 1:num_vars
        eval([active_vars{j} '_new=[' active_vars{j} '_new;' active_vars{j} '];']);
    end
end

%Copy *_new back to the original variable names
for i = 1:num_vars
    eval([active_vars{i} '=' active_vars{i} '_new;']);
end
    
%Save the concatenated variables to the specified location
save(output, active_vars{:})

end