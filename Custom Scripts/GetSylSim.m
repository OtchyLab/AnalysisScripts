function GetSylSim(File, DOB)
%Takes in a single .mat file containing the variables:
%   FeatureSet
%   FeatureKeys
%   FeatureTotals
%   SylTypes

%Pre=allocate know vars
SimStats = [];

%Load the specified file
load(File);

%Determine date of recording and bird age to save with output
RecDate = regexp(File, '_', 'split'); %Parse file name for date
RecDateNum = datenum(RecDate(2), 'yyyymmdd');
DOBDateNum = datenum(DOB, 'yyyymmdd');

%Determine the unique syllables ID's in the files
UniqueSyls = sort(unique(SylTypes(SylTypes>0 & SylTypes<100)));
numSyls = length(UniqueSyls);

for i = 1:numSyls
    FeatureSubSet = FeatureSet(SylTypes==UniqueSyls(i),:);
    FeatureSubSet = ScaleFeats(FeatureSubSet);
    SimStats(i).FeatureKeysSubSet = FeatureKeys(SylTypes==UniqueSyls(i));
    [SimStats(i).Short_sq, SimStats(i).Long_sq] = PSucker(FeatureSubSet, FeatureSubSet);
    SimStats(i).SylTyp = UniqueSyls(i);
    SimStats(i).Short = 1-squareform(SimStats(i).Short_sq);
    SimStats(i).Long = 1-squareform(SimStats(i).Long_sq);
    SimStats(i).Short_m = mean(SimStats(i).Short);
    SimStats(i).Long_m = mean(SimStats(i).Long);
    SimStats(i).RecDate = char(RecDate(2));
    SimStats(i).BirdAge = RecDateNum-DOBDateNum;
end

%Save SimStats to file
filename = [char(RecDate(1)) '_' char(RecDate(2)) '_SylSim.mat'];
save(filename, 'SimStats');

end