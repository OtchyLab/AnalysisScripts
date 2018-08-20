
%Establish vars
similarity = [];
harmonicIdx = [];
files = [];

%Setup ages and ranges to process
binAges = 50:10:120;
bounds = [45, 56;...
                57, 63;...
                67, 73;...
                77, 83;...
                87, 93;...
                97, 103;...
                107, 113;...
                117, 123];

%Cycle through all matching files in the folder
filenames = dir('*SimScore*.mat');
idx = 1;
%This loop points to files in the directory
for i = 1:length(filenames)
    %Load out data
    load(filenames(i).name)

    %This loop points to syllables within a single file
    for j = 1:size(meandata,1)
        
        %This loop points to intervals of time defined in bounds
        for k = 1:size(bounds,1)
            %Create a logical index in which true means the timepoint is within the age range
            mask = Agedata>bounds(k,1) & Agedata<bounds(k,2);
            
            %Retrieve any similarity scores that were in that range
            picks = meandata(j,mask);
            
            %Capture similarity score
            if isempty(picks)
                %If there's nothing in there, mark it as NaN
                similarity(idx,k) = NaN;
            else
                %Otherwise, 
                similarity(idx,k) = picks;
            end
        end
        
        %Things we want to keep
    harmonicIdx(idx) = harmonic(j);
    files{idx} = filenames(i).name;
    
    %Update index
    idx = idx + 1;
    end


    
    
    
    
    
    
end