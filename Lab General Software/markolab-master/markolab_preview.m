function AX=markolab_preview(DATA,DATA_FS,AUDIO,AUDIO_FS)
%
%
%

order=.1;
spect_size=1;
max_fs=8e3;

AX=[];

AX(1)=subplot(spect_size+1,1,1:spect_size);
[s,f,t]=zftftb_pretty_sonogram(AUDIO,AUDIO_FS,'filtering',300,'norm_amp',1,'clipping',[-3 2]);
imagesc(t,f/1e3,s);axis xy;
ylim([0 8]);
colormap(jet);
ylabel('Fs (kHz)');
ylimits=ylim();
set(AX(1),'YTick',ylimits,'XTick',[],'TickDir','out','TickLength',[0 0]);

timevec=[1:size(DATA,1)]/DATA_FS;

AX(2)=subplot(spect_size+1,1,spect_size+1);
plot(timevec,DATA);
linkaxes(AX,'x');
xlabel('Time (s)');
ylimits=ylim();
set(AX(2),'TickDir','out','TickLength',[0 0],'YTick',ylimits);
set(zoom(gcf),'ActionPostCallback',{@preview_change_zoom,AX,order});
set(pan(gcf),'ActionPostCallback',{@preview_change_zoom,AX,order});

end

function preview_change_zoom(gcbo,eventdata,handles,order)

for i=1:length(handles)
	ylimits=ylim(handles(i));
	%ylimits=round(ylimits/order)*order;
	%ylim(handles,[ylimits]);
	set(handles(i),'TickDir','out','TickLength',[ 0 0 ],'YTick',ylimits);
end

end
