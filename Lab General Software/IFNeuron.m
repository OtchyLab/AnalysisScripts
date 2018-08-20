%Adapting Integrate and fire from Christoff Koch chapter 14
%voltages in millivolts, currents in nanoamps, 
% resistances in megaohms
% capacitance in nanoFarads
%
%fully vectorized version


clear all
%Number of neurons
nNeuron = 10;
%timing
dt = .05; %millisecond time step
tEnd = 150; %maximum simulation time
%build the time vector
time = 0:dt:tEnd ;


%membrane parameters
tMemb = 3*ones(1,nNeuron);  %mSec time constant
rMemb = 40*ones(1,nNeuron); %memb resistance
cMemb = tMemb/rMemb *ones(1,nNeuron); %memb capacitance
vRest = 0;  %millivolts resting potential
tRef = 2.5*ones(1,nNeuron); %mSec: length of the absolute refactory period. v is held at zero for this time

vThresh = 16*ones(1,nNeuron); %threshold voltage

%adapting (K) cond: gAdapt
tAdapt = 20*ones(1,nNeuron); %mSec time constant of slow recovery after an AP 
Ginc = .0/dt*ones(1,nNeuron);  %.05 microSiemans: Instant increase in gAdapt when a spike occurs 

%synaptic connections
%weight matrix --
%each column represents the input weights to one neuron
%e.g. element(1,2) is the weight connecting neuron 1 to neuron 2

weight = zeros(nNeuron,nNeuron);

%weight = 3 * randn(nNeuron,nNeuron) + .0;
%weight = 2 * rand(nNeuron,nNeuron) - .25;

%synaptic current rises instantly, then decays exponentially
%decay rate of synaptic current 
Irate = .3*ones(1,nNeuron);
Idecay = 1 - Irate*dt;

%init state variables
vMemb = zeros(1,nNeuron);
spike = zeros(1,nNeuron);
gAdapt = zeros(1,nNeuron);
spikeTime = -inf*ones(1,nNeuron); %time of occurence of last spike

%initial input currents
Io = 0.25*ones(1,nNeuron);
r = 1.4;
Is = 3*[1 1 1 1 1 r r r r r]; %steady current
Ir = zeros(1,nNeuron);
%step the time
tindex = 1; %pointer for output arrays

%init plotting variables
Vout = zeros(length(time),nNeuron) ;  
Aout = zeros(length(time),nNeuron) ;
Sout = zeros(length(time),nNeuron) ;


for t=time

    %input current to cells
    %Istate =  spike*weight + Istate .* Idecay ; %decaying synaptic currents
    Ir = 0.01*randn(size(Is)) + .99*Ir ;
    I =  Is .* (t>50 & t<100).*(1 + 1*Ir) + Io.*(1 + 2.0*Ir) ;
    
    %add to plotting variables
    Vout(tindex,:) = vMemb;  
    Aout(tindex,:) = I ;
    Sout(tindex,:) = spike;

    
    %update state variables
    %use second order Runge-Kutta
    % dgAdapt/dt = (Ginc.*spike - gAdapt)./tAdapt
    % dv/dt = 

    k1g = dt*((Ginc.*spike - gAdapt./tAdapt)) ;

    k1v = dt*( I./cMemb - vMemb.*(1 + rMemb.*gAdapt)./tMemb) .* ...
        ((t-spikeTime) >= tRef) ;

    k2g = dt*(Ginc.*spike - (gAdapt + 0.5*k1g)./tAdapt) ;
    k2v = dt*( I./cMemb - (vMemb + 0.5*k1v).*(1 + rMemb.*(gAdapt + 0.5*k1g))./tMemb) .* ...
        ((t+0.5*dt-spikeTime) >= tRef) ;

   
    gAdapt = gAdapt + k2g;
    vMemb = vMemb + k2v;

    %generate the spike
    spike = (vMemb >= vThresh);
    %and its time of occurence
    spikeTime(spike) = t;
    %hit the membrane potential
    vMemb(spike) = 0.;

    %update time index
    Tout(tindex) = t;
    tindex = tindex+1;

end %main time step for-loop

figure(1)
clf
np=3;

subplot(np,1,1)
plot(Tout,Vout)
ylabel('V')

subplot(np,1,2)
plot(Tout,Aout)
ylabel('Iin')

subplot(np,1,3)
plot(Tout,Sout)
ylabel('Spikes')
set(gca,'ylim',[-.2 1.5])
xlabel('time (mS)')

spiketrain = Sout(fix(50/dt):fix(100/dt),:);

figure(2)
clf
pn=1;
% subplot(pn,1,1)
% plot(spiketrain)
% xlabel('Time (samples)');

%compute times of spikes in seconds into a CELL ARRAY
for i=1:nNeuron
    tSpike{i} = .001*time(find(spiketrain(:,i)>0.5))+.05 ;
end

subplot(pn,1,1)

for i=1:nNeuron
    for j=i+1:nNeuron
        
        col = [1 0 0]*(i<=5)*(j<=5)+[0 0 1]*(i>5)*(j>5);
        
        %Compute the distances at different time scales
        for s=1:13
            TimeScale(s) = .001 * 1.5^(s-1);
            
            %calc J.Victor distance
            %cost = 1/TimeScale(s) ; %per Sec
            %dVictor(s) = spkd(tSpike{i},tSpike{j},cost) ;
            
            %calc vanRossum dist
            %make dt depend on timescale
            tc = TimeScale(s) ; % Sec
            dt = tc/100 ;
            tk = 0:dt:4*tc;
            %convert tend to seconds
            t = 0:dt:tEnd*.001;
            
            spiketrain1 = zeros(size(t));
            spiketrain2 = zeros(size(t));
            spiketrain1(fix(tSpike{i}/dt)) = 1;
            spiketrain2(fix(tSpike{j}/dt)) = 1;
            
            ekernel = exp(-tk/tc);
            cSpike1 = conv(ekernel,spiketrain1);
            cSpike2 = conv(ekernel,spiketrain2);
            dRossum(s) =  2*dt/tc * sum((cSpike1-cSpike2).^2) ;
        end
        
        
        %semilogx(TimeScale,dVictor)
        
        semilogx(TimeScale,dRossum,'color',col) 
        hold on
        drawnow
        
    end
end


xlabel('Time scale (sec)')
ylabel('Distance')
%legend('Victor Dist','Rossum Dist')
