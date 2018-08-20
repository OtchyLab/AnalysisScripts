function identifySyllablesInExper(exper, filenums)

for(filenum = filenums)
    loadAudio
    call identifySyllables
    
    save syllable information:
        datastructure?
            syllable 
            identity
            dir
            filename
            starttime_absolute 
            endtime_absolute
            startndx in audio file
            endndx in audio file
            prevent overlap
            exper
            filenum
            
end
            

