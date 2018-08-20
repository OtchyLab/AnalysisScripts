function changeDrugStatus(exper, start, finish)

[file,path] = uigetfile('*.mat', 'Choose an audio annotation .mat file to load:');
if isequal(file,0) || isequal(path,0)
else
    annotationFileName = [path,file];
    audioAnnotation = aaLoadHashtable(annotationFileName);
    

end  
for filenum=start:finish
filename = getExperDatafile(exper,filenum,0);
    if(audioAnnotation.containsKey(filename))
            currAnnot = audioAnnotation.get(filename);
            currAnnot.drugstatus='Muscimol in LMAN';
            currAnnot.drugindex=2;
            audioAnnotation.put(filename, currAnnot);
    end

end
aaSaveHashtable(annotationFileName, audioAnnotation);
