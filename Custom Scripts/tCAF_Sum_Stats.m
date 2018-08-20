function [output] = tCAF_Sum_Stats()
%Function to output the stats for the Ali et al 2013 papaer (3/20/13)

%Request and load first day dataset
[output.fname1,path1] = uigetfile('C:\Users\Tim\Desktop\Ali Again\*.mat');
load([path1 output.fname1],'intervals')
output.day1_intervals = intervals;
clear('intervals')

%Request and load second day dataset
[output.fname2,path2] = uigetfile('C:\Users\Tim\Desktop\Ali Again\*.mat');
load([path2 output.fname2],'intervals')
output.day2_intervals = intervals;
clear('intervals')

%Pur626
%output.target_ints = [5,6];
%output.nontarget_ints = [1,2;3,4;7,7];

%Pur692
% output.target_ints = [5,6];
% output.nontarget_ints = [1,2;3,4;7,7];

%Pur755
output.target_ints = [3,4];
output.nontarget_ints = [1,2;5,6;7,7];
output.days = 3;

int_diff = output.day2_intervals.m-output.day1_intervals.m;

output.tCAFoverall = sum(int_diff(output.target_ints))/output.days; %tCAF overall per day
output.tCAFsylonly = int_diff(output.target_ints(1))/output.days; %tCAF syl only change per day
output.tCAFgaponly = int_diff(output.target_ints(2))/output.days; %tCAF gap only change per day

output.tCAFgapStart = output.day1_intervals.m(output.target_ints(2)); %tCAF gap starting duration
output.tCAFgapEnd = output.day2_intervals.m(output.target_ints(2)); %tCAF gap ending duration

output.tCAFsylStart = output.day1_intervals.m(output.target_ints(1)); %tCAF syl starting duration
output.tCAFsylEnd = output.day2_intervals.m(output.target_ints(1)); %tCAF syl ending duration

output.tCAFnonTarget = mean(sum(int_diff(output.nontarget_ints),2))/output.days; %tCAF nontarget overall per day

[sname,path] = uiputfile('C:\Users\Tim\Desktop\Ali Again\*.mat');
save([path, sname],'output');