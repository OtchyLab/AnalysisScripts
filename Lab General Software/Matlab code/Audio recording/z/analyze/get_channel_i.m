function y=get_microhpone_i(yAInfo,name)
y=-1;
for i=1:length(yAInfo.ObjInfo.Channel)
    if strcmp(yAInfo.ObjInfo.Channel(i).ChannelName,name)
      y=i;
  end    
 end