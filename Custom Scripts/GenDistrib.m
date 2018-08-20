function [Output] = GenDistrib(Input)
%Takes in a Feature File from GetWaveFeatures and produces an output
%structure with ecdf and ecdf arrays for AM, FM, Entropy, and Pitch.  The
%Output structure has the following organization:
%   Output.epdf
%       Output.epdf.AMx
%       Output.epdf.AMy
%       Output.epdf.FMx
%       Output.epdf.FMy
%       Output.epdf.Entx
%       Output.epdf.Enty
%       Output.epdf.GPx
%       Output.epdf.GPy
%   Output.ecdf
%       Output.ecdf.AMx
%       Output.ecdf.AMy
%       Output.ecdf.FMx
%       Output.ecdf.FMy
%       Output.ecdf.Entx
%       Output.ecdf.Enty
%       Output.ecdf.GPx
%       Output.ecdf.GPy

%Define Output structure
%Output.epdf
      Output.epdf.AMx = [];
      Output.epdf.AMy = [];
      Output.epdf.FMx = [];
      Output.epdf.FMy = [];
      Output.epdf.Entx = [];
      Output.epdf.Enty = [];
      Output.epdf.GPx = [];
      Output.epdf.GPy = [];
%Output.ecdf
      Output.ecdf.AMx = [];
      Output.ecdf.AMy = [];
      Output.ecdf.FMx = [];
      Output.ecdf.FMy = [];
      Output.ecdf.Entx = [];
      Output.ecdf.Enty = [];
      Output.ecdf.GPx = [];
      Output.ecdf.GPy = [];
%Output.ecdf.yproj
      Output.ecdf.yproj.AM = [];
      Output.ecdf.yproj.FM = [];
      Output.ecdf.yproj.Ent = [];
      Output.ecdf.yproj.GP = [];

%option to use select feature file with 
if nargin == 0
    [file, folder]=uigetfile('*.mat', 'Feature File');
    Input = importdata([folder file]);
end

%Grab all epdf arrays
[Output.epdf.AMx,Output.epdf.AMy]=epdfplot(Input(1,:),100);
[Output.epdf.FMx,Output.epdf.FMy]=epdfplot(Input(2,:),100);
[Output.epdf.Entx,Output.epdf.Enty]=epdfplot(Input(3,:),100);
[Output.epdf.GPx,Output.epdf.GPy]=epdfplot(Input(5,:),100);

%Grab all ecdf arrays
[Output.ecdf.AMy,Output.ecdf.AMx]=ecdf(Input(1,:));
[Output.ecdf.FMy,Output.ecdf.FMx]=ecdf(Input(2,:));
[Output.ecdf.Enty,Output.ecdf.Entx]=ecdf(Input(3,:));
[Output.ecdf.GPy,Output.ecdf.GPx]=ecdf(Input(5,:));

%Project ecdf onto the y axis at 100 evenly spaced points
steps = 0:.01:1;
array = [];

for i=1:length(steps)
    index = find(Output.ecdf.AMy>=steps(i),1,'first');
    Output.ecdf.yproj.AM(i) = Output.ecdf.AMx(index);
    
    index = find(Output.ecdf.FMy>=steps(i),1,'first');
    Output.ecdf.yproj.FM(i) = Output.ecdf.FMx(index);
    
    index = find(Output.ecdf.Enty>=steps(i),1,'first');
    Output.ecdf.yproj.Ent(i) = Output.ecdf.Entx(index);
    
    index = find(Output.ecdf.GPy>=steps(i),1,'first');
    Output.ecdf.yproj.GP(i) = Output.ecdf.GPx(index);
end


end