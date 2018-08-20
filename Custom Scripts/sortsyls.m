function sortsyls
%Prompt user for directory
dname = uigetdir('C:\');

cd(dname); %change to directory)
filelist=dir('*.wav'); %get file list

for i=1:length(filelist)
   fname = filelist(i).name;
   parse = regexp(fname, '_', 'split');
   syltype=parse(end);
   syltype = regexprep(syltype, '.wav', '');
   if ~isdir([dname '\' char(syltype)])
       mkdir([dname '\' char(syltype)])
   end
   movefile(fname,[dname '\' char(syltype) '\' fname])
end

end