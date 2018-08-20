% this might not have to be it's own function ... p simple
function name = formFileName(birdname, SOMETHING, Y, M, D, H, MN, S)
    name = strcat(birdname, "_", SOMETHING, "_", num2str(Y), "_", num2str(M), ...
    "_", num2str(D), "_", num2str(H), "_", num2str(MN), "_", num2str(S)); 
end