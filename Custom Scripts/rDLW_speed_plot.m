%Creates figure showing a comparison of rDLW resolution and speed to prior
%work

clear

%Data from the Saha paper (all in log units)
pnts = [-6.10812897676737,2.08499408677294;-4.83370637162433,2.17026184355849;-4.73401229560271,2.35424628346349;-4.42779509935199,2.34681809740698;-4.15138838661754,2.67653181902593;-2.19005404982749,2.90791981468631;-2.11215583552433,2.80069981331796;-2.00894314505488,2.70882488051372;-1.74934759023780,2.35653338285457;-1.75345264042693,2.67836931768202;-1.91081289767674,3.01541348606727;-0.252763578430898,3.02434685719312;2.47308748643867,3.31762337141908];

%rDLW data (mm3/hr, nm res)
rDLW = [1.6089, 2.845];

spd = pnts(:,1);
res = pnts(:,2);

%Plot the figure
figure (33); clf

old = loglog(10.^spd, (10.^res)*2.2, '.k', 'MarkerSize', 14); hold on
new = loglog(10.^rDLW(1), 10.^rDLW(2), 'sqg', 'MarkerSize', 14);


%format
set(gca, 'Box', 'off', 'TickDir', 'out')

xlim([1e-7,1e3])
ylim([1e2,1e4])

xlabel('Volumetric processing rate (mm3/hr)')
ylabel('Resolution (nm)')








