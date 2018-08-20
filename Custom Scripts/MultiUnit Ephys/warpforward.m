function [mwarpout,swarpout] = warpforward(master, slave, fs)

%Master pattern to warp everything to
warper = master{1};

if length(master)==length(slave)
    for i=1:length(master)
        mwarpee=master{i}; %remove cell structure
        swarpee=slave{i}; %remove cell structure
        [mwarpout{i},p,q]=dtw(warper,mwarpee,fs); %warp master timeseries to pattern
        swarpout{i} = warpit(warper,swarpee,fs,p,q); %warp slave timeseries to previous
    end
else
    print('Master and slave datasets are not the same length. Exiting...')
end



end