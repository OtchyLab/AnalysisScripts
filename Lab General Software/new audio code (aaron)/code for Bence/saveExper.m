function saveExper(exper)
%The exper contains a directory folder specifying where the files are
%saved.  This direction needs to be updated if the files are moved.
%Rootdir is the directory in which the bird folder resides.

save([exper.dir,'exper.mat'],'exper'); 