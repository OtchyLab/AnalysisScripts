function truncateAnnotation(from,to,newAnnot, keys, elements)
for i=from:to
    newAnnot.put(keys{i},elements{i})
end