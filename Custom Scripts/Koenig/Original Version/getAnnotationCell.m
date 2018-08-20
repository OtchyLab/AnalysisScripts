function fieldVector=getAnnotationCell(handles,field)

for k = 1:size(handles,2)
    AnnStruct=handles{1,k};
    fieldVector{k} = getfield(AnnStruct, field);
end