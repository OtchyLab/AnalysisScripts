%%%%%% draw the stuff 
samples_available=get(ai,'SamplesAvailable');
if samples_available>0
    y=peekdata(ai,samples_available)';
    %%%%%% compress y
    ly=size(y,2);   
    zf=floor(ly/points_per_pixel/screenSize(3));y=y(:,1:zf:end);
    %y=interp1(max(y,0),1:zf:ly,'linear')
    %z=interp1(min(y,0),1:zf:ly,'linear');
    %ly=size(y,2);
    for i=1:length(ai.Channel)
        
        set(hLine(i),'xdata',(1:zf:ly)/scanrate,'ydata',y(i,:));
        
        set(hAxes(i),'XLim',[0 ly/scanrate],'YLim',[min(y(i,:)) max(y(i,:))]);
        
        if micro>-1
            set(hTrig,'xdata',soft_trig_time*[1 1],'ydata',[min(y(1,:)) max(y(1,:))]);
        end   
    end
    set(h_rep_fun,'xdata',samp_times,'ydata',-ones(1,length(samp_times)));
   samp_times;
else 
    %etime(clock,ai.InitialTriggerTime);
    fprintf('samples_available=0\n');
end
