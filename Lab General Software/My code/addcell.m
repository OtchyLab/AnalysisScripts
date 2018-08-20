function totalCell=addcell(totalCell,presentCell, startrec, endrec,cellno)

pointer=length(totalCell.rec);
totalCell.rec(pointer+1:pointer+(endrec-startrec)+1)=cellno;
totalCell.motif(pointer+1:pointer+(endrec-startrec)+1,1:2)=presentCell.motif(startrec:endrec,1:2);
totalCell.syll(pointer+1:pointer+(endrec-startrec)+1,1:size(presentCell.syll,2))=presentCell.syll(startrec:endrec,1:size(presentCell.syll,2));
totalCell.spikes(pointer+1:pointer+(endrec-startrec)+1,1:size(presentCell.spikes,2))=presentCell.spikes(startrec:endrec,1:size(presentCell.spikes,2));

