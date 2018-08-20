
close_analyze_window; %delete any windows that might be there
threshold=0.85;
online_default=0;
do_plot_defaultA=1;
%songtype_prefixA_default='D'; % 'D'= Directed, 'U'='Undirected', 'B'='BOS'
loadpath_default='c:\Documents and Settings\Emberbence\My Documents\Temporary Data Files\Matlab files\';
callerA_default='Feedback';

define_dependent_variablesA;
loadpath=loadpath_default;
%%%%%%%%%%%%%%%%%%%%% GET BIRDNAMES
birdnames=struct2cell(dir(loadpath));
if size(birdnames,2)<3 
    fprintf('No Data\n');
    break
end
birdnames=birdnames(:,3:end);
[birdname_dates birdname_dates_i]=sort(birdnames(2,:)); %sorts the dates the files were modifies
birdnames=birdnames(1,birdname_dates_i); %the oldest file is first
birdnameA=birdnames{end};
[dirnames,dirnameA]=get_directories(loadpath,birdnameA); %dirnames is the list of dates, dirameA is the last date
get_filenames;

control_analyze_figure;
%%%%%%%%%%%%%%%%%%%%% GET FILES and display them
%directed_vs_undirectedA;

if do_plotA
    control_analyze_figure_plot;
   
if length(namehelpA)>0
 analyze_plot;
end
end
