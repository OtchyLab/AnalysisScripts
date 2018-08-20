% New Pipeline for eventual FinchScope --> talonPipeline conversion

% 1) Run FS_AV_Parse to generate .MAT files for all of the videos
% 2) Run deleteFrames to get rid of the first 10 frames of the video
%    (changes the .MAT file and the .CSV file)
% 3) Run FS_DFF_STD_Image
% 4) Run FS_BatchDff
% 5) Take an ROI cast from 