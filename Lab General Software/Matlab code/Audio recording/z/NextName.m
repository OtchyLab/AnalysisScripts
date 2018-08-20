function y=NameHelp(namehelp)
zerosh={'000','00','0',''};
i=str2num(namehelp(end-3:end)); i=i+1;
t=num2str(i);
t=[zerosh{length(t)} t];
y=[namehelp(1:end-4) t];