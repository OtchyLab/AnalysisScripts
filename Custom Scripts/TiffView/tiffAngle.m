% % 7/18/18
% % Script that will hopefully find the angle of tilt for the tif files...
% % that was a bad description IDK how to say it in a smart way yet

% % Step 1: Get file list
% % set directory 
%%
directory = '/Users/sadiela/Desktop/TestCubeTiffs';
% get names of all tiff files in directory 
files = dir(fullfile(directory, '*.tif'));
[names, index] = sortrows({files.name}');

% Step 2: Read in series of images 
tiff.im = [];

% this reads in all of the slices of a file
for j = 1:length(names)
    i = 0;
    try
        while true
            i = i + 1;
            tiff(j).im(i,:,:) = imread(strcat(directory, "/", names{j}), i); % dis might be d issue: handles.files is just a cell array of strings
            %imagesc(squeeze(tiff(i,:,:)));
            % when does this loop end?
        end
    catch ME
    end
end

filename = 'alltiffs.mat';
save(filename, 'tiff');

%% SHOULDN'T HAVE TO RUN ANYTHING ABOVE HERE, JUST LOAD FILE


% load in all tiff images 
load('/Users/sadiela/Desktop/alltiffs.mat');

% select tiff to work on: TestCube_V-1.06_00001_00001.tif
curr_tiff = tiff(7).im;

% the image is 512 pixels, but its 666 microns across, so each pixel is
% 1.301 microns 
% so, at some point, multiply all of the x&y values by 1.301 in order to 
% get them to the same scale as the z values! 


% select only images that are transitioning from dark to light and cut off
% the edges
rel_curr_tiff_edgeless = curr_tiff(40:52, 11:502, 11:502);
inspect = squeeze(rel_curr_tiff_edgeless(4,:,:));
imagesc(inspect);

% convert to binary image 
% for now, take cuttoff be 5500 (adjust later)
for i = 1:13
    for j = 1:492
        for k = 1:492
            if(rel_curr_tiff_edgeless(i,j,k) <= 5500)
                rel_curr_tiff_edgeless(i,j,k) = 0;
            else
                rel_curr_tiff_edgeless(i,j,k) = 1;
            end
        end
    end
end

% Write a function to find the range of files to use for angle calculation!
% (based on # of light and dark pixels

% look at first and last image and see total:

first = squeeze(rel_curr_tiff_edgeless(1,:,:));
last = squeeze(rel_curr_tiff_edgeless(13,:,:));

%%
% I BROKE SOMETHING :((( TIFF ISN'T RIGHT ANYMOREEEE


% these will be filled with all the points along the light/dark boundary
% for each image
allx = [];
ally = [];
allz = [];

% loop through each image to get the x & y points 
for j = 1:13
    % second end goal: struct with x & y REGRESSION vals for each image
    curImage = squeeze(rel_curr_tiff_edgeless(j,:,:)); % get current image
    [rows, cols] = size(curImage); 
    [startRow, startCol] = findBottomPix(curImage); % this wont work in every case!
    points = zeros(2, 1); % this will hold the x&y coordinates for the first "white" pixel in each row
    counter = 1;
    for i = 1:rows
        col = startCol; % this will mess everything up if the black slants left
        while (curImage(i, col) == 1 && col < cols)
            col = col+1;
        end
        if (col ~= cols) 
            points(1,counter) = i; % y components (row #)
            points(2,counter) = col; % x components (column #)
            counter = counter+1;
        %else
            %points(1,i) = 0;
            %points(2,i) = 0;
        end
    end
    x = fliplr(points(2,:));
    y = points(1,:);
    if(x(length(x)) ~= 1 && length(x) ~=1)
        % dont call linReg if corrcoef doesn't work... annoying
        % switch to microns now?
        x = x * 1.301; % is this the right place to do this? It is likely the problem with the final angle value
        y = y * 1.301;
        [xreg, yreg] = simpleLinReg(x,y);
        allx = [allx xreg];
        ally = [ally yreg];
        zadd = zeros(1, length(xreg)) + j;
        allz = [allz zadd];
    end
   
    [r, m, b] = regression(x, y);
    regressStats(1,j) = r;
    regressStats(2,j) = m;
    regressStats(3,j) = b;
end

allz = allz(1:length(allx));

%%

% plot all points
scatter3(allx, ally, allz);

% perform planar regression
targPlane = planarRegression(allx,ally, allz);
A = plane(1);
B = plane(2);
C = plane(3);

% Ax + By + Cz + D = 0
D = dot([allx(1) ally(1) allz(1)],[A B C]);

targPlane = [A B C];


[U1, V1] = angleAdjustments('/Users/sadiela/Desktop/TestCubeTiffs', 'TestCube_baseline_00001_00001.tif')
[U2, V2] = angleAdjustments('/Users/sadiela/Desktop/TestCubeTiffs', 'TestCube_U-1.64_V2.49_00001_00001.tif')
[U3, V3] = angleAdjustments('/Users/sadiela/Desktop/TestCubeTiffs', 'TestCube_U-2.16_V-1.11_00001_00001.tif')
[U4, V4] = angleAdjustments('/Users/sadiela/Desktop/TestCubeTiffs', 'TestCube_U-3.46_00001_00001.tif')
[U5, V5] = angleAdjustments('/Users/sadiela/Desktop/TestCubeTiffs', 'TestCube_U2.59_V-0.79_00001_00001.tif')
[U6, V6] = angleAdjustments('/Users/sadiela/Desktop/TestCubeTiffs', 'TestCube_U4.58_00001_00001.tif.tif')
[U7, V7] = angleAdjustments('/Users/sadiela/Desktop/TestCubeTiffs', 'TestCube_V-1.06_00001_00001')
[U8, V8] = angleAdjustments('/Users/sadiela/Desktop/TestCubeTiffs', 'TestCube_V3.02_00001_00001.tif')






xyPlane = [0 0 1]; % normal vectors
xzPlane = [0 1 0];
yzPlane = [1 0 0];

%%%%%%%%%%%%%%%%%%%%%%

% V direction: (left right)
Vline = lineBetweenPlanes(targPlane, xzPlane);
xaxis = [1 0 0];
Vangle = angleBetweenVectors(Vline, xaxis);


% U direction: (up down)
Uline = lineBetweenPlanes(targPlane, yzPlane);
yaxis = [0 1 0];
Uangle = angleBetweenVectors(Uline, yaxis);

% one should be almost zero and it is so yay


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% More garbage 

% create x, y
[xx,yy] = ndgrid(1:10, 1:10);
% calculate corresponding z

z = (-A*xx - B*yy -D)/C;

surf(xx,yy,z);

hold;

% linspace x and y points
xmax = max(allx);
xmin = min(allx);
% calculate z based on equation
% plot with mesh?
z = (A*x + B*y)/C;

patch(xpoints, ypoints, zpoints, [1 1 0]);


% XY plane is [0 0 1 0]
targPlane = [A B C];

% find the angle between the two planes!
ang = angleBetweenPlanes(targPlane, xyPlane);
% convert to degrees
Uang_d = (360 * Uangle)/(2*pi); % 179.651 degrees (.3490 degrees)
Vang_d = (360 * Vangle)/(2*pi); % 

% the target angle for this specific image was 1.06 degrees in the V
% direction, so something is wrong

%% 

% If I want to decompose the angle into its u and v components, I need to
% find its projection onto the xz and yz planes (see drawings)
% these will correspond to the u and v directions


%%

% SCRAPS 

dark = squeeze(rel_curr_tiff_edgeless(1,:,:));
high_dark_val = max(max(dark));

middle = squeeze(rel_curr_tiff_edgeless(5,:,:));


%white = squeeze(rel_curr_tiff_edgeless(20,:,:));
min_light_val = min(min(white));

imshow(dark)
imshow(middle)
imshow(white)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rows, cols] = size(middle);
[startRow, startCol] = findBottomPix(middle);
points = zeros(2, 1); % this will hold the x&y coordinates for the first "white' pixel in each row
counter = 1;
for i = 1:rows
    col = startCol;
    while (middle(i, col) == 1 && col < cols)
        col = col+1;
    end
    if (col ~= cols) 
        points(1,counter) = i; % y components (row #)
        points(2,counter) = col; % x components (column #)
        counter = counter+1;
    %else
        %points(1,i) = 0;
        %points(2,i) = 0;
    end
end


x = fliplr(points(2,:));
y = points(1,:);


xbar = mean(x);
xstd = std(x);
ybar = mean(y);
ystd = std(y);
rmat = corrcoef(x,y); % I hate matlab
r = 0.9831;

b = r * (ystd/xstd);
a = ybar - b*xbar;

% Y = bX + a
% 

X_reg = 1:length(x);

Y_reg = b*X_reg + a;


%%
% SCRAPS
for j = 1:492
    for k = 1:492
        if(rightTiff(j,k) <= 5500)
            rightTiff(j,k) = 0;
        else
            rightTiff(j,k) = 1;
        end
    end
end