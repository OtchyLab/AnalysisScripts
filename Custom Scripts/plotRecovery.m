%Plot the dataOut structure

%Load from file
outputName =  'C:\Users\Tim\Desktop\Temp Nif Data\RecoveryData.mat';
% load(outputName,'dataOut')

numSets = length(dataOut);

h(1) = figure(1);
cla; hold on

h(2) = figure(2);
cla; hold on

h(3) = figure(3);
cla; hold on

for i = 1:numSets
    n = [dataOut(i).birdName, ' ', dataOut(i).birdChan];
    if strcmp(dataOut(i).birdName, 'Grn046')
        col = 'r';
    elseif strcmp(dataOut(i).birdName, 'Grn121')
        col = 'b';
    elseif strcmp(dataOut(i).birdName, 'Grn141')
        col = 'g';
    end
    figure(1)
    errorbar(dataOut(i).days, dataOut(i).recoveryPowNorm.m, dataOut(i).recoveryPowNorm.std, ['-' col 'x'], 'DisplayName', n)
    title('Power Recovery')
    
    figure(2)
    errorbar(dataOut(i).days, dataOut(i).recoveryDistNorm.m, dataOut(i).recoveryDistNorm.std, ['-' col 'x'], 'DisplayName', n)
    title('Distance Recovery')
    
    figure(3)
    errorbar(dataOut(i).days, dataOut(i).recoveryCorrNorm.m, dataOut(i).recoveryCorrNorm.std, ['-' col 'x'], 'DisplayName', n)
    title('Correlation Recovery')
end

%Lesioned birds
bird1PowM =mean([dataOut(1).recoveryPowNorm.m; dataOut(2).recoveryPowNorm.m; dataOut(3).recoveryPowNorm.m]); %Grn121
bird1PowS =std([dataOut(1).recoveryPowNorm.m; dataOut(2).recoveryPowNorm.m; dataOut(3).recoveryPowNorm.m]);
bird1Days = dataOut(1).days;

bird2PowM =mean([dataOut(4).recoveryPowNorm.m; dataOut(5).recoveryPowNorm.m; dataOut(6).recoveryPowNorm.m]); %Grn141
bird2PowS =std([dataOut(4).recoveryPowNorm.m; dataOut(5).recoveryPowNorm.m; dataOut(6).recoveryPowNorm.m]);
bird2Days = dataOut(4).days;

NaN(2, max([bird1Days, bird2Days]))

%Intact birds
bird4PowM = [dataOut(7).recoveryPowNorm.m]; %Pur683
bird4PowS =[dataOut(7).recoveryPowNorm.std];
bird4Days = dataOut(7).days;

bird5PowM =mean([dataOut(8).recoveryPowNorm.m; dataOut(9).recoveryPowNorm.m; dataOut(10).recoveryPowNorm.m]); %Pur690
bird5PowS =std([dataOut(8).recoveryPowNorm.m; dataOut(9).recoveryPowNorm.m; dataOut(10).recoveryPowNorm.m]);
bird5Days = dataOut(8).days;

bird6PowM =mean([dataOut(11).recoveryPowNorm.m; dataOut(12).recoveryPowNorm.m]); %Pur696
bird6PowS =std([dataOut(11).recoveryPowNorm.m; dataOut(12).recoveryPowNorm.m]);
bird6Days = dataOut(11).days;

bird7PowM =mean([dataOut(13).recoveryPowNorm.m; dataOut(14).recoveryPowNorm.m]); %Pur696
bird7PowS =std([dataOut(13).recoveryPowNorm.m; dataOut(14).recoveryPowNorm.m]);
bird7Days = dataOut(13).days;



