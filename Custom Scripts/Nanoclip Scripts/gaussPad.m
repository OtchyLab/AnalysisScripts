function padmat = gaussPad(mat)
%Take in a matrix and pad it with gaussians to help define time
% Assume that columns define the time dimension
%

% Written by TMO; last mod 09/12/2018

bDebug = false;

%Gaussian kernal definition
width = 11; %gaussian width in units of the input matrix (for specs, this is typically ms)
kern = gausswin(width);

%Size of input matrix
[~, n] = size(mat);

%Create impulse matrix marking number and timing of gaussians
step = 9; %gaussian spacing in units of the input matrix (for specs, this is typically ms)
wing = floor(step/2);
numGs = floor(n/step);
Imat = zeros(numGs, n);

for i = 1:numGs
    idx = 1 + (step*(i-1));
    Imat(i, idx) = 1;
    
    q = conv(kern, Imat(i,:));
    Imat(i,:) = q(wing:n+wing-1);
    
end

%Copy to output var
padmat = Imat;

%Debugging check
if bDebug
    figure; imagesc(Imat);
end
