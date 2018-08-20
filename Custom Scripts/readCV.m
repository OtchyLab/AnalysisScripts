function readCV
%Function reads the text file from produiced by teh Gamry Potentistat in
%cyclic voltemmetry mode.
%
%Code taken pretty much entirely from Yarden's github code:
%https://github.com/gardner-lab/lab-operations/blob/master/utilities/Electrochemistry%20measurements/ReadCV.m
%
%Last modified by TMO (02/06/2017)

% Pad dimensions for CIC calculations
pW = 0.007; %in cm
pH = 0.0035; %in cm
pArea = pH*pW; %in cm^2

%Ask user which file to load; check for suffix
[filename, pathname, filterindex] = uigetfile('X:\CV*.DTA', 'Pick a DTA file');

%Read DTA script
a = DTAread(fullfile(pathname,filename),'\t',0);

%Parse out the sections we care about
title = a{1};
a = a(2:end);

figure(100); clf
legnd={};
clrs=parula(numel(a)-1);
smth = 5;
for cycnum=2:numel(a)
    plot(a{cycnum}.data(:,3),smooth(a{cycnum}.data(:,4),smth),'Color',clrs(cycnum-1,:));
    hold on;
    legnd{cycnum-1}=['cyc #' num2str(cycnum-1)];
end
xlim([-1, 1]); ylim([-.2e-6, .2e-6]);
set(gca,'FontSize',14);
ylabel('I (amp)');
xlabel('V vs Ag|AgCl');
legend(legnd);

%CIC calculation on the last trace
x = a{cycnum-1}.data(:,3);
y = smooth(a{cycnum-1}.data(:,4),smth);
y(y>0) = 0;
CIC = polyarea(x,y*1000/pArea);
t = ['CIC = ' num2str(CIC,4) 'mC/cm^2'];
text(-0.75, 1e-7, t, 'FontSize', 14)
