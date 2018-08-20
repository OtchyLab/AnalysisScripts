function exper = createExper(rootdir)
%Creates a folder for all files related to this experiment.  Also saves a
%.mat file to this folder containing the experiment description.

if(~exist('rootdir'))
    rootdir = pwd;
end

birdname = input('Enter a bird name: (no spaces or strange characters)','s');
birddesc = input('Enter a description of the bird:','s');
expername = input('Enter a experiment name:','s');
experdesc = input('Enter a description of the exper:','s');

mkdir(rootdir, birdname);
mkdir([rootdir,'/',birdname], expername);

exper.dir = [rootdir,'\',birdname,'\',expername,'\'];
exper.birdname = birdname;
exper.birddesc = birddesc;
exper.expername = expername;
exper.experdesc = experdesc;
exper.datecreated = datestr(now,30);

exper.desiredInSampRate = input('Enter the desired input sampling rate:');
exper.audioCh = input('What hw channel will audio be on: (-1 if no audio)');
exper.sigCh = input('Enter vector of other hw channels to be recorded: ([] if none)');

save([rootdir,'\',birdname,'\',expername,'\exper.mat'], 'exper');
