function simulateRA(varargin)

% we ignore the voltage dependence of the currents since reversal potential
% is far from subthreshold potentials.
if ~isempty(varargin)
    HVC=varargin{1};
    if (nargin==2)
        age=varargin{2};
    end
end
time_step=1; %ms
time_sim=960; %length of simulation in ms
meanHVCsyn=0.1;
trials=40;
bin_size=3; %3 ms for psth bins
sdHVCsyn=0.3;
fraction=0.95; %ca. fraction of active HVC neurons
V_rev_AMPA=0;
V_rev_NMDA=0;
gainHVC=0.05;
gainLMAN=0.04;
gainNoise=10;
Vrest=-75;
V_AMPA=0;
V_NMDA=0;
fraction_AMPA_in_LMAN=0.15;
fraction_AMPA_in_HVC=0.6;
cMg=1;%mg concentration in mM
Vreset=-60;
gainInh=3.5;
LMAN_fr=50;
Vthresh=-50;
noHVC=floor(time_sim)/12; %each neuron fires for 12 ms, every 3 ms

tau_AMPA=20/time_step;
tau_NMDA=100/time_step;
for i=1:(time_sim/time_step)
    filter_AMPA(i)=exp(-(i-1)/tau_AMPA);
    filter_NMDA(i)=exp(-(i-1)/tau_NMDA);
end
if (~exist('HVC'))
    HVC=zeros(noHVC,floor(time_sim/time_step));
    W_HVC = random('Normal',meanHVCsyn,sdHVCsyn,1,noHVC); %weights for the different HVC neurons

    on_off=rand(1,noHVC);
    W_HVC(find(on_off>=fraction))=0;

    W_HVC(find(W_HVC<0))=0;
    W_HVC(find(W_HVC>1))=1;

    for i=1:noHVC %spiking of HVC neurons
        HVC(i,(i-1)*(12/time_step)+1:(3/time_step):i*(12/time_step))=W_HVC(1,i);
    end
end
if (exist('age'))
    for i=1:noHVC %spiking of HVC neurons
        if (HVC(i,(i-1)*(12/time_step)+1)>0.37)
            HVC(i,(i-1)*(12/time_step)+1:(3/time_step):i*(12/time_step))=HVC(i,(i-1)*(12/time_step)+1)+age*(1-HVC(i,(i-1)*(12/time_step)+1));
        else
            HVC(i,(i-1)*(12/time_step)+1:(3/time_step):i*12)=HVC(i,(i-1)*(12/time_step)+1)-age*(HVC(i,(i-1)*(12/time_step)+1));
        end
    end
end
% figure;plot(LMAN);
figure;imagesc(HVC);

for i=2:floor(time_sim/time_step)   
    [row,col]=find(HVC(:,1:i-1)>0);
    for j=1:length(row)
        HVC_cond(j,i) = gainHVC*(HVC(row(j),col(j))*(fraction_AMPA_in_HVC*filter_AMPA(i-col(j))+(1-fraction_AMPA_in_HVC)*filter_AMPA(i-col(j)))); 
    end
end

HVC_sumCond=sum(HVC_cond,1);
figure;
ax(1)=subplot('Position',[0.06 0.4 0.9 0.55]);
%subplot(2,1,1);
for t=1:2*trials
    V(1)=Vrest;
    noise=rand(1,floor(time_sim/time_step));
    noise=noise-0.5;
    for i=1:floor(time_sim/time_step)
        if(rand<(0.001*LMAN_fr))
            LMAN(i)=1;
        else
            LMAN(i)=0;
        end
    end
    
    LMAN_cond_NMDA=[];
    LMAN_cond_AMPA=[];
    for i=2:floor(time_sim/time_step)   
        LMAN_spikes=find(LMAN(1:i-1)>0);
        for j=1:length(LMAN_spikes)
            LMAN_cond_NMDA(j,i) = gainLMAN*(filter_NMDA(i-LMAN_spikes(j))); 
            LMAN_cond_AMPA(j,i) = gainLMAN*fraction_AMPA_in_LMAN*(filter_AMPA(i-LMAN_spikes(j)));
        end
    end

    LMAN_sumCond_NMDA=sum(LMAN_cond_NMDA,1);
    LMAN_sumCond_AMPA=sum(LMAN_cond_AMPA,1);
    if(t<trials+1)
        LMAN_on=1;
        Inh=gainInh*(gainHVC+gainLMAN*(2.5*fraction_AMPA_in_LMAN+0.25*(1-fraction_AMPA_in_LMAN)));
    else
        LMAN_on=0;
        Inh=gainInh*(gainHVC);
    end
    spike(1:floor(time_sim/time_step))=0;
    for i=2:floor(time_sim/time_step) 
        NMDA_V=1/(1+exp(-0.062*V(i-1))*(cMg/3.57));
        I_HVC(i)=HVC_sumCond(i)*(V(i-1)-V_AMPA);
        I_inh(i)=(V(i-1)-Vrest)*Inh;
        I_LMAN_NMDA(i)=LMAN_sumCond_NMDA(i)*(V(i-1)-V_NMDA)*NMDA_V*LMAN_on;
        I_LMAN_AMPA(i)=LMAN_sumCond_AMPA(i)*(V(i-1)-V_AMPA)*LMAN_on;
        if (V(i-1)-I_HVC(i)-I_inh(i)-I_LMAN_NMDA(i)-I_LMAN_AMPA(i)>Vthresh)
            V(i)=Vreset;
            spike(i)=1;
        else
            V(i)=V(i-1)-I_HVC(i)-I_inh(i)-I_LMAN_NMDA(i)-I_LMAN_AMPA(i)+noise(i);
        end

    end
    spikes=(find(spike>0)/(1000/time_step)); %spike time in seconds
    for l=1:length(spikes)
        line([spikes(l) spikes(l)], [t-1 t-0.2]);
    end
    allSpikes{t}=spikes;
end
line([0 time_sim/(1000*time_step)],[trials trials],'color','r');
drugSpikes=cat(2,allSpikes{1:trials});
noDrugSpikes=cat(2,allSpikes{trials+1:2*trials});
for i=1:floor(time_sim/bin_size)
    bin(i)=(i*bin_size)/1000;
end
psth_drug = histc(drugSpikes,bin)/(trials*(bin_size/1000));
psth_noDrug=histc(noDrugSpikes,bin)/(trials*(bin_size/1000));
hold off;
ax(2)=subplot('Position',[0.06 0.06 0.9 0.28]);
linkaxes(ax,'x');
plot(bin,psth_noDrug);
hold on;
plot(bin,psth_drug,'color','r');
hold off
% figure;plot(V);
% hold on;
% plot(I_HVC,'color','k');
% plot(I_LMAN_NMDA,'color','c');
% plot(I_LMAN_AMPA,'color','r');
% plot(I_inh,'color','g');
% hold off;
button = questdlg('save HVC?');
if (strcmp(button,'Yes'))
    [FileName,PathName]=uiputfile('*.mat','Save HVC As:','HVC');
    save([PathName FileName],'HVC');
end

% I=I_HVC+I_LMAN+I_Leakage;
% deltaV=deltaT*I/C;
