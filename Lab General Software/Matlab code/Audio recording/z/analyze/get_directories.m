function [y,z]=get_directories(loadpath,birdnameA);
%%%%%%%%%%%%%%%%%%%%% GET DIRECTORIES --- DATES 
dirnames=struct2cell(dir([loadpath birdnameA '/20*']));
if size(dirnames,2)==0
    y='';
else
    y=sort(dirnames(1,1:end));
end
if length(y)>0
   z=y{end};
else
    z='';
end
