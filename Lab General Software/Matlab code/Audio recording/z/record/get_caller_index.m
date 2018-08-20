function y=get_caller_index(callerfield,name)
y=-1;
for i=1:length(callerfield)
    if strcmp(callerfield{i},name)
        y=i;
    end
end
