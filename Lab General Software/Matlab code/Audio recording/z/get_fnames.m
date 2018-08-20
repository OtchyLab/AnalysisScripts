function y=get_fnames(savepath,birdname,dirname,name_prefix,songtype_prefix,caller_prefix,bytes) 
d=dir([savepath birdname '\' dirname '\' name_prefix birdname songtype_prefix '-' caller_prefix '*']);
ci=0; fnames={};fhelp={};
for i=1:size(d,1)
    if d(i).bytes>bytes
%%%%%%% the following is a solution to the problem that the dir function
%%%%%%% does not distinguish betwen upper case and lower case
          if length(caller_prefix)==0 | strcmp(caller_prefix,d(i).name(end-7-length(caller_prefix):end-8))
            ci=ci+1; fnames{ci}=d(i).name; fhelp{ci}=d(i).name(end-7:end-4);
        end    
    end
end
[yh i]=sort(fhelp); %fhelp is simply the numbers
y=fnames(i); %what we get is the filenames in order