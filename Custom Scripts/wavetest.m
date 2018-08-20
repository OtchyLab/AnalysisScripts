function wavetest

%Misc Constants (for now...)
threshold=25;
min_stop=20/1000;
min_dur=10/1000;
buffer=5/1000;
direct = uigetdir(pwd,'Directory name');

%Batch Calls
files = dir('*.wav');
for k = 1:numel(files);
   file = files(k).name;

%Just check the in and out
[handles.wav handles.fs]=wavread(file);

%Plot waveform
subplot(3,1,1);
plot((0:length(handles.wav)-1)/handles.fs,handles.wav,'c');
yl = max(abs(handles.wav))*1.05;
ylim([-yl yl]);
ylabel('Amplitude (ADU)')
xlim([0 (length(handles.wav)-1)/handles.fs]);
set(gca,'color',[0 0 0],'xtick',[],'ytick',[]);
hold off

% Plot sound amplitude
handles.amplitude = 1000*abs(handles.wav);
subplot(3,1,2);
plot((0:length(handles.wav)-1)/handles.fs,handles.amplitude,'b');
hold on
handles.thresline = plot([0 (length(handles.wav)-1)/handles.fs],repmat(threshold,1,2),':r');
hold off
box off
ylim([0 max(handles.amplitude)+eps]);
ylabel('Amplitude (ADU)')
xlim([0 (length(handles.wav)-1)/handles.fs]);

%Segment into syllables
% Find threshold crossing points
handles.segments = [];
a = [0; handles.amplitude; 0];
handles.segments(:,1) = find(a(1:end-1)<threshold & a(2:end)>=threshold)-1;
handles.segments(:,2) = find(a(1:end-1)>=threshold & a(2:end)<threshold)-1;
a = a(2:end-1);

% Eliminate short intervals
i = [find(handles.segments(2:end,1)-handles.segments(1:end-1,2) > min_stop*handles.fs); length(handles.segments)];
handles.segments = [handles.segments([1; i(1:end-1)+1],1) handles.segments(i,2)];

% Eliminate short syllables
i = find(handles.segments(:,2)-handles.segments(:,1)>min_dur*handles.fs);
handles.segments = handles.segments(i,:);

% Attach buffer
handles.segments(:,1) = handles.segments(:,1) - buffer*handles.fs;
handles.segments(find(handles.segments(:,1)<1),1) = 1;
handles.segments(:,2) = handles.segments(:,2) + buffer*handles.fs;
handles.segments(find(handles.segments(:,1)>length(a)),2) = length(a);

%Plot segment overlay
handles.selectedseg = ones(1,size(handles.segments,1));
handles.segmenttitles = cell(1,size(handles.segments,1));
    for c = 1:length(handles.segmenttitles)
        handles.segmenttitles{c} = '';
    end
subplot(3,1,3)
cla
set(gca,'color',get(gcf,'color'),'xcolor',get(gcf,'color'),'ycolor',get(gcf,'color'));
handles.segmentplots = zeros(1,size(handles.segments,1));
handles.segmentlabels = zeros(1,size(handles.segments,1));
for c = 1:size(handles.segments,1)
    plot([handles.segments(c,1) handles.segments(c,2)]/handles.fs,[0 0],'linewidth',4);
    hold on;
    % plot([handles.segments(c,1) handles.segments(c,2)]/handles.fs,[0 0],cols{handles.selectedseg(c)+1},'linewidth',4);
    handles.segmentlabels(c) = text((handles.segments(c,1)+handles.segments(c,2))/handles.fs/2, .75, handles.segmenttitles{c});
end
ylim([-1 1]);
xlim([0 (length(handles.wav)-1)/handles.fs]);
hold off

%Save all syllables as individual wavs
seg_list = find(handles.selectedseg==1);
isin =size(handles.segments,1);
%direct = uigetdir(pwd,'Directory name');
for c = seg_list
    fname = [file '_syll_' num2str(c,'%03.f') '.wav'];
    wavwrite(handles.wav(round(handles.segments(c,1):min([handles.segments(c,2) length(handles.wav)]))),handles.fs,16,[direct '\' fname]);
end
hold off
box off
%msgbox([num2str(length(seg_list)) ' syllables written.'],'Done');

end