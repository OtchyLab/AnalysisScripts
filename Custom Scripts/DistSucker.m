function [Shorty, Longy] = DistSucker(List1, List2)
%Function will take in the lists of syllable features in List1 and List2,
%and calulcate the distance between each row in List1 and every row in
%List2.  The features in the input vars should be (Mx33)and have already
%been scaled by ScaleFeats.

%Test input array for formatting
[m, n] = size(List1);
[o, p] = size(List2);

if n ~=33
    error('List1 not appropriately formatted.  Must be Mx33.')
end

if p ~=33
    error('List2 not appropriately formatted.  Must be Ox33.')
end

%Pre-allocate the needed variables
Shorty = zeros(m,o,4,4);
Longy = zeros(m,o);
Temp =[];

%Call HowFar to populate the Distance Matrices
for i=1:m
    for j=1:o
        [Dshort, Longy(i,j)] = HowFar(List1(i,:), List2(j,:));
        %Shorty(i,j,:,:) = Dshort;
        Temp = [Temp;Dshort(1,:)';Dshort(2,:)';Dshort(3,:)';Dshort(4,:)'];
    end
end
Shorty = Temp;
end