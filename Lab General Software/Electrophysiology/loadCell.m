function cell = loadCell(cellNum,cellDir)

if(~exist('cellDir'))
    cellDir = 'C:\MATLAB6p5\work\Aaron\Electrophysiology\Cells\';
end

filename = getCellFilename(cellNum);

load([cellDir,filename]);
