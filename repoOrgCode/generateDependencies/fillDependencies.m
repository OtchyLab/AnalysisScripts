%%% This script takes a list of folder names that correspond with different
%%% steps of analysis (e.g. SongBlaster, TweetVision) that is assigned to
%%% the variable folderNames. It generates a unique (nonrepeating) list of
%%% dependent files for all of the functional scripts in each of these
%%% folders and puts these files into the target directory. 
% Note that it is assumed that all folders for which you are generating
% dependencies are in the baseDirectory, defined below. It is also assumed
% that the "functional" scripts in the folders are those with the same
% name (and subsequent iterations of the names, such as SapSimilarity2.m
% and SapSimilarity3.m)
% This script calls the functions: getPathNames, getRelScripts, and
% generateDependenceList
baseDirectory = "/Users/labperson/Desktop/Software/Custom Scripts";
folderNames = {"Songblaster","TweetVision","TweetVisionLite","StretchEm", ...
    "StretchEmLite","Metermeter","AnnotDisplay","FinchShocker","Koenig",...
    "SapSimilarity","Talon","Wave_clus"};

% generate pathnames for each folder
folderPaths = getPathNames(folderNames, baseDirectory,1);
namesPaths = getPathNames(folderNames, baseDirectory,2);

% generate list of relevant files
filenames = getRelScripts(namesPaths);

% generate list of dependencies and put them in a "Common file" folder
target = '/Users/labperson/Desktop/Software/Common Files';
generateDependenceList(filenames, folderPaths, target);

