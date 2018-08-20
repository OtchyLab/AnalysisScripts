function [Features] = GetWaveFeatures(song,fs)
%Function takes a .wav file as input and returns the spectral features
%calculated by Sound Analysis for Matlab.  Features is an mxn matrix where
%n are timepoints (typically covering ~1ms of sound) and m are the various
%spectral features:
%   
%   1: m_AM
%   2: m_FM
%   3: m_Entropy
%   4: m_amplitude
%   5: m_PitchGoodness
%   6: m_Pitch
% 

%option to leave out the default sampling rate
if nargin == 1
    fs = 44100;
end

%option to use default sampling rate AND select song file with GUI
if nargin == 0
    fs = 44100;
    [file, folder]=uigetfile('*.wav', 'Song File');
    song = [folder file];
end

%Define output matrix
Features = [];

%read in single from file
wav = wavread(song);


%execute function to calculate spectral components
[m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude ,gravity_center, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight ]=deriv(wav,fs);

%Collapse relevant data to single output variable
Features = [m_AM'; m_FM'; m_Entropy'; m_amplitude'; m_PitchGoodness'; m_Pitch'];






end