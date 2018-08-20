%files is the output of the 'dir' function, and is a list of all files in a
%directory. These files may be in random temporal order. The function sorts the 
%files according to their date and
%returns the index for the file (file_no) that was chosen (in terms of temporal 
%order within the directory). The file_order
%provides the index to the original list of files in the proper order. 

function [file_order,file_no]=sortfiles(files, filename) 

date=zeros(length(files)-2,length(files(1).date));
for i=1:length(files)-2; %there are always two elements devoted to 'dots', thus subtract 2.
    date(i,:)=files(i+2).date;
    if (strcmp(filename, files(i+2).name))
        file_index=i;
    end
end
[spacer,file_order]=sortrows(date);%file_order = pointer to the files in the right temporal order
file_order=file_order+2; %because of the two 'dot'elements
[spacer,file_sort]=sort(file_order);
file_no=file_sort(file_index)+2; %file_no: the chosen file in terms of temporal order

%i=index