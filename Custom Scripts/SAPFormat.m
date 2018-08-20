function [Bird, DataOut, VarOuts, NetVar] = SAPFormat(DOB)
%Function is designed to format the output matrices of the SAP batch
%similarity function.  Importing that output matrix produces two MATLAB
%variables -- one text, one numerical data.  In the text variable, the
%column of interest should have the format:
%       BirdNum_dRecnum_YYYYMMDDTHHMMSSchan#_syll_###_SylType.wav
%of which, BirdNum, the date, and SylType are the most important.  The
%output of this function is a string vector with the bird's name and a
%numerical vector with:
%       Bird's age, the syllable type, and the similarity score

%Read in DOB as a MM-DD-YYYY and calculate the datenum
DOBnum = datenum(DOB, 'mm-dd-yyyy');

%Get the .xls data output from SAP batch program
[FileName,PathName] = uigetfile('*.xlsx');
cd(PathName);
[Data, Text] = xlsread(FileName, 'Sheet1');

%For the DAT-stripped file structure
Text=Text(2:end,3); %Strip down to a single useful column (Song2 Title)
Data = Data(:,4); %Strip down to a single useful column (Accuracy)

%For the WAV file structure:
% Text=Text(:,3); %Strip down to a single useful column (Song2 Title)
% Data = Data(:,2); %Strip down to a single useful column (Accuracy)


%Filter out all entries for which Accuracy=100% (same file comparison) or
%Accuracy<=0 (file read empty)
TextFiltered = Text(Data<100 & Data>0);
DataFiltered = Data(Data<100 & Data>0);

%Create variables
Bird = [];
DataOut = [];
VarOuts = []; %Age, SylType, Var, StD
NetVar = [];
Letters = ['ABCDEFGHIJK'];
letters = ['abcdefghijk'];

%Break out the text data and extract relevant values
for i=1:length(TextFiltered)
%For the DAT-stripped file structure 
    Parts=regexp(char(TextFiltered(i)), '_', 'split');
    Bird{i} = Parts(1); %Bird Name
    Datestr = char(Parts(3));
    Datestr = Datestr(1:8);
    DataOut(i,1) = round(datenum(Datestr, 'yyyymmdd') - DOBnum); %Age
    SylStr = char(Parts(6));
    DataOut(i,2) = str2num(SylStr(1)); %SylType
    DataOut(i,3) = DataFiltered(i); %Similarity

% %For the WAV file structure
%     Parts=regexp(char(TextFiltered(i)), '_', 'split');
%     Bird{i} = Parts(1); %Bird Name
%     Datestr = [char(Parts(3)) char(Parts(4)) char(Parts(5))];
%     %Datestr = Datestr(1:8);
%     DataOut(i,1) = round(datenum(Datestr, 'yyyymmdd') - DOBnum); %Age
%     SylStr = char(Parts(11));
%     SylStr = SylStr(1);
%     index = strfind(Letters,SylStr); %SylType
%     if isempty(index)
%         index = strfind(letters,SylStr);
%     end
%     DataOut(i,2) = index;
%     DataOut(i,3) = DataFiltered(i); %Similarity
end
Bird = Bird'; %Transpose for ease of reading

%Grab all of the unique values from output vectors
uAges = sort(unique(DataOut(:,1)));
uSyls = sort(unique(DataOut(:,2)));

%For every unique Age and Syllable, calculate variability and Std
m=1;
for i=1:length(uAges)
    SubSet = DataOut(DataOut(:,1)==uAges(i),:);
    n=m;
    for j=1:length(uSyls)
        SubSet2 = SubSet(SubSet(:,2)==uSyls(j),:);
        if ~isempty(SubSet2)
            VarOuts(m,1) = uAges(i); %Age
            VarOuts(m,2) = uSyls(j); %SylType
            VarOuts(m,3) = mean(SubSet2(:,3)); %mean Similarity
            VarOuts(m,4) = ((100-VarOuts(m,3))/50)*100; %Var
            VarOuts(m,5) = std(SubSet2(:,3)); %StD
            m=m+1;
        end
    end
    NetVar = [NetVar; uAges(i), mean(VarOuts(n:m-1,3)), std(VarOuts(n:m-1,3))];
end


end

