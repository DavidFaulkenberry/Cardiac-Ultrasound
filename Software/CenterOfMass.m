%% BME543 Final Project
% Group 4: Tyler Meathe, David Faulkenberry, Connor Johnson
% Detecting Center of Mass in Left Ventricle
% Assisted by Cooper Moore and Olaf Von Ramm
 
%% Read in Data
M = readDicom3D('patient2.dcm')
framenumber=M.NumVolumes;
 
% Extract actual dimensions in cm of the image
w_cm=M.widthspan;
h_cm=M.heightspan;
d_cm=M.depthspan;
 
% # of samples per dimension in the original data
h=M.height;
w=M.width;
d=M.depth;
 
% Distance between each sample point in each dimension
h_dis=h_cm/h;
w_dis=w_cm/w;
d_dis=d_cm/d;
density=1.05; % g/mL density of blood

%% Filter that Andy Used
% Filter each frame's volume in the height dimension (1D convolution)
for frame = 1:framenumber
    M.dataHeightFiltered(:,:,:,frame) = FilterHeightwise(M.data(:,:,:,frame));
end
 
%%  Create 3D matrix with respect to individual time frame
for time=1:framenumber
%Tested using data from the first time frame
%V=squeeze(M.data(:,:,:,time)); %Frame raw data
V=squeeze(M.dataHeightFiltered(:,:,:,time)); %Frame raw data
V=double(V); %Convert it to double
 
ns= 208; %number of samples in each dimension of the new cube data matrix
 
%% Resampling Data into New Matrix
rmax= w;
rpmax=ns;
cmax= h;
cpmax= ns;
dmax= d;
dpmax=ns;
NewHeart= zeros(ns,ns,ns); % Create the new cube data matrix in 3D
 
 
for i=1:ns
    r=((i-1).*(rmax-1)./(rpmax-1))+1;
    rref=floor(r);
    rtop=ceil(r);
    a= r-rref;
    for j=1:ns
        c=((j-1).*(cmax-1)./(cpmax-1))+1;
        cref=floor(c);
        ctop=ceil(c);
        b= c-cref;
        for k=1:ns
            Bright= V(rref,cref,k).*(1-a).*(1-b)...
                +V(rref,ctop,k).*a.*(1-b)+...
                V(rtop,cref,k).*(1-a).*(b)+...
                V(rtop,ctop,k).*a.*b;
            NewHeart(i,j,k)= Bright;
        end
    end
end
 
%% Viewing Original slices
Normalize=NewHeart./max(max(NewHeart));
xyview = squeeze(Normalize(:,:,ns./2));
 
xzview=squeeze(Normalize(:,ns./2,:));
xzview=imrotate(xzview,270);
 
figure(1);
imshowpair(xyview,xzview,'montage')
title('Middle Slice in Horizontal and Vertical Directions, Unfiltered');
xlabel(['Frame ',num2str(time)])
 
%% Choose Long Apex and Base Points
figure(2);
imshow(xzview)
title('Select Apex then Base of LV for Long Axis Measurement');
xlabel(['Frame ',num2str(time)])
[xz, yz] = ginput(2);
LongAxisTop = [xz(1), yz(1)];
LongAxisBottom = [xz(2), yz(2)];

longaxisdist=longaxisdistance(LongAxisTop(1),LongAxisTop(2),LongAxisBottom(1),LongAxisBottom(2));
 
Storage(time).longaxis=longaxisdist;
 
%% Slice Creation
slice4 = round((yz(2) + yz(1))/2);
slice2 = round((slice4 + yz(1))/2);
slice1 = round((yz(1) + slice2)/2);
slice3 = round((slice2 + slice4)/2);
slice6 = round((yz(2) + slice4)/2);
slice5 = round((slice4 + slice6)/2);
slice8 = round(yz(2));
slice7 = round((slice8 + slice6)/2);
 
slice=[slice1 slice2 slice3 slice4 slice5 slice6 slice7 slice8];
Storage(time).slice=slice;
 
%% slice number in the depth plane
for s=1:length(slice)
% Averaging 5 Slices in xy and xz
thickness = 8;
i = 0;
ThiccsliceXY = 0;
ThiccsliceXZ = 0;
 
%% Slice 1
if s==1
for i=0:thickness
    NewSliceXY = NewHeart(:,:,slice(s) + i);
    ThiccsliceXY = ThiccsliceXY + NewSliceXY;
end
ThiccsliceXY = ThiccsliceXY./ thickness;
 
 
% normalizing averaged slices
ThiccsliceXY=ThiccsliceXY./max(max(ThiccsliceXY));
Filtered2Dxy = zeros(ns);
Filtered2Dxz = zeros(ns);
 
for x = 1:ns
    for y = 1:ns
        %filter xy
        if ThiccsliceXY(x, y) <= 0.67
            Filtered2Dxy(x, y) = 0;
        else
            Filtered2Dxy(x, y) = ThiccsliceXY(x, y);
        end
    end
end
 
figure(7)
imshow(ThiccsliceXY)
title('Thick Slice Before Filtering');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])
 
%% Displaying Thresholded Slices
figure(3)
imshow(Filtered2Dxy);
title('Slice Filtered by Thresholding (Removes Minimum Brightness)');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])
 
%% Filtering 
figure(4);
Ixy=imbinarize(Filtered2Dxy,'global');
edgexy=edge(Ixy,'log');
imshow(Ixy);
title('Filtering Threshold Slice with Otsu Approximation');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])

figure(5)
imshow(edgexy)
title('Edge Detection from Otsu Approximation');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])
 
%% Crop Images to Remove Noise
edgexy = Crop(edgexy, 90, 140, 80, 140, ns);
 
%% Find Circle
figure(6);
imshow(edgexy)
title('Edge Detection to Generate Circles - Image Has Been Cropped');
[centers,radii] = imfindcircles(edgexy,[8 250]) ;
indradii=find(max(radii)==radii);
for i=1:length(indradii)
    newradii(i)=radii(indradii(i));
    for j=1:2
        newcenters(i,j)= centers(indradii(i),j);
    end
end

%% Area of Circle with 8 Pie Slices (Every 45 Degrees)
DegreeRadius = eightradii(zeros(1, 8), newcenters, edgexy, ns, newradii);
Storage(time).newradii.slice1=DegreeRadius;
 
%% Starting Storage
Storage(time).circlecenter.slice1=newcenters;
viscircles(newcenters, newradii,'Color','b');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])
newradii=newradii.*(w_cm./ns);
Storage(time).circleradii.slice1=newradii;
CircleArea=pi.*(newradii).^2;
CircleVolume=CircleArea.*longaxisdist
Storage(time).circlearea.slice1=CircleArea;
Storage(time).circlevol.slice1=CircleVolume;
Storage(time).zcenter.slice1=(yz(1)+slice1)/2;
 
%% Area of a Slice
Area = 0;
nonzeros = find(DegreeRadius ~= 0);
count = length(nonzeros);
for i=1:count
    Area = Area + (1/count).*pi.*(DegreeRadius(i).^2);
end
 
Storage(time).newarea.slice1=Area;
Volume = Area.*longaxisdist./count;
MassSlice=density.*Volume;
Storage(time).newvolume.slice1=Volume;
Storage(time).mass.slice1=MassSlice;


end
 
%% Slice 2
if s==2
for i=0:thickness
    NewSliceXY = NewHeart(:,:,slice(s)+ i);
    ThiccsliceXY = ThiccsliceXY + NewSliceXY;
end
ThiccsliceXY = ThiccsliceXY./ thickness;
 
 
% normalizing averaged slices
ThiccsliceXY=ThiccsliceXY./max(max(ThiccsliceXY));
Filtered2Dxy = zeros(ns);
Filtered2Dxz = zeros(ns);
 
for x = 1:ns
    for y = 1:ns
        %filter xy
        if ThiccsliceXY(x, y) <= 0.4650
            Filtered2Dxy(x, y) = 0;
        else
            Filtered2Dxy(x, y) = ThiccsliceXY(x, y);
        end
    end
end
 
%% Displaying Thresholded Slices
figure(3)
imshow(Filtered2Dxy);
title('Slice Filtered by Thresholding (Removes Minimum Brightness)');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])
 
%% Filtering 
figure(4);
Ixy=imbinarize(Filtered2Dxy,'global');
edgexy=edge(Ixy,'log');
imshow(Ixy);
title('Filtering Threshold Slice with Otsu Approximation');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])

figure(5)
imshow(edgexy)
title('Edge Detection from Otsu Approximation');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])
 
%% Crop Images to Remove Noise
edgexy = Crop(edgexy, 83, 140, 70, 140, ns);
 
%% Find Circle
figure(6);
imshow(edgexy)
title('Edge Detection to Generate Circles - Image Has Been Cropped');
[centers,radii] = imfindcircles(edgexy,[8 250]) ;
indradii=find(max(radii)==radii);
for i=1:length(indradii)
    newradii(i)=radii(indradii(i));
    for j=1:2
        newcenters(i,j)= centers(indradii(i),j);
    end
end
 
%% Area of Circle with 8 Pie Slices (Every 45 Degrees)
DegreeRadius = eightradii(zeros(1, 8), newcenters, edgexy, ns, newradii);
Storage(time).newradii.slice2=DegreeRadius;
 
%% Starting Storage
Storage(time).circlecenter.slice2=newcenters;
viscircles(newcenters, newradii,'Color','b');
title(['Frame ',num2str(time),', Slice ', num2str(s)])
newradii=newradii.*(w_cm./ns);
Storage(time).circleradii.slice2=newradii;
CircleArea=pi.*(newradii).^2;
CircleVolume=CircleArea.*longaxisdist
Storage(time).circlearea.slice2=CircleArea;
Storage(time).circlevol.slice2=CircleVolume;
Storage(time).zcenter.slice2=(slice2+slice1)/2;
 
%% Area of a Slice
Area = 0;
nonzeros = find(DegreeRadius ~= 0);
count = length(nonzeros);
for i=1:count
    Area = Area + (1/count).*pi.*(DegreeRadius(i).^2);
end
 
Storage(time).newarea.slice2=Area;
Volume = Area.*longaxisdist./count;
MassSlice=density.*Volume;
Storage(time).newvolume.slice2=Volume;
Storage(time).mass.slice2=MassSlice;
end
 
%% Slice 3
if s==3
for i=0:thickness
    NewSliceXY = NewHeart(:,:,slice(s) + i);
    ThiccsliceXY = ThiccsliceXY + NewSliceXY;
end
ThiccsliceXY = ThiccsliceXY./ thickness;
 
 
% normalizing averaged slices
ThiccsliceXY=ThiccsliceXY./max(max(ThiccsliceXY));
Filtered2Dxy = zeros(ns);
Filtered2Dxz = zeros(ns);
 
for x = 1:ns
    for y = 1:ns
        %filter xy
        if ThiccsliceXY(x, y) <= 0.36
            Filtered2Dxy(x, y) = 0;
        else
            Filtered2Dxy(x, y) = ThiccsliceXY(x, y);
        end
    end
end
 
figure(7)
imshow(ThiccsliceXY)
title('Thick Slice Before Filtering');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])
 
%% Displaying Thresholded Slices
figure(3)
imshow(Filtered2Dxy);
title('Slice Filtered by Thresholding (Removes Minimum Brightness)');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])
 
%% Filtering 
figure(4);
Ixy=imbinarize(Filtered2Dxy,'global');
edgexy=edge(Ixy,'log');
imshow(Ixy);
title('Filtering Threshold Slice with Otsu Approximation');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])

figure(5)
imshow(edgexy)
title('Edge Detection from Otsu Approximation');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])
 
%% Crop Images to Remove Noise
edgexy = Crop(edgexy, 80, 140, 74, 140, ns); %top bottom left right
 
%% Find Circle
figure(6);
imshow(edgexy)
title('Edge Detection to Generate Circles - Image Has Been Cropped');
[centers,radii] = imfindcircles(edgexy,[8 250]) ;
indradii=find(max(radii)==radii);
for i=1:length(indradii)
    newradii(i)=radii(indradii(i));
    for j=1:2
        newcenters(i,j)= centers(indradii(i),j);
    end
end
%% Area of Circle with 8 Pie Slices (Every 45 Degrees)
DegreeRadius = eightradii(zeros(1, 8), newcenters, edgexy, ns, newradii);
Storage(time).newradii.slice3=DegreeRadius;
 
%% Starting Storage
Storage(time).circlecenter.slice3=newcenters;
viscircles(newcenters, newradii,'Color','b');
title(['Frame ',num2str(time),', Slice ', num2str(s)])
newradii=newradii.*(w_cm./ns);
Storage(time).circleradii.slice3=newradii;
CircleArea=pi.*(newradii).^2;
CircleVolume=CircleArea.*longaxisdist
Storage(time).circlearea.slice3=CircleArea;
Storage(time).circlevol.slice3=CircleVolume;
Storage(time).zcenter.slice3=(slice2+slice3)/2;
%% Area of a Slice
Area = 0;
nonzeros = find(DegreeRadius ~= 0);
count = length(nonzeros);
for i=1:count
    Area = Area + (1/count).*pi.*(DegreeRadius(i).^2);
end
 
Storage(time).newarea.slice3=Area;
Volume = Area.*longaxisdist./count;
MassSlice=density.*Volume;
Storage(time).newvolume.slice3=Volume;
Storage(time).mass.slice3=MassSlice;
end
 
%% Slice 4
if s==4
for i=0:thickness
    NewSliceXY = NewHeart(:,:,slice(s) + i);
    ThiccsliceXY = ThiccsliceXY + NewSliceXY;
end
ThiccsliceXY = ThiccsliceXY./ thickness;
 
% normalizing averaged slices
ThiccsliceXY=ThiccsliceXY./max(max(ThiccsliceXY));
 
%% Attempting to Threshold Filter in 2D
Filtered2Dxy = zeros(ns);
Filtered2Dxz = zeros(ns);
 
for x = 1:ns
    for y = 1:ns
        %filter xy
        if ThiccsliceXY(x, y) <= 0.30
            Filtered2Dxy(x, y) = 0;
        else
            Filtered2Dxy(x, y) = ThiccsliceXY(x, y);
        end
    end
end
 
%% Displaying Thresholded Slices
figure(7)
Filtered2Dxz = imrotate(Filtered2Dxz,270);
imshowpair(Filtered2Dxy,Filtered2Dxz,'montage');
 
%% Filtering 
Ixy=imbinarize(Filtered2Dxy,'global');
edgexy=edge(Ixy,'Roberts');
imshowpair(Ixy,edgexy,'montage');
 
%% Crop Images to remove noise
edgexy = Crop(edgexy, 85, 145, 65, 145, ns);
 
%% Find Circle
figure(6);
imshow(edgexy)
title('Edge Detection to Generate Circles - Image Has Been Cropped');
[centers,radii] = imfindcircles(edgexy,[8 250]) ;
indradii=find(max(radii)==radii);
for i=1:length(indradii)
    newradii(i)=radii(indradii(i));
    for j=1:2
        newcenters(i,j)= centers(indradii(i),j);
    end
end
 
%% Area of Circle with 8 Pie Slices (Every 45 Degrees)
DegreeRadius = eightradii(zeros(1, 8), newcenters, edgexy, ns, newradii);
Storage(time).newradii.slice4=DegreeRadius;
%% Starting Storage
Storage(time).circlecenter.slice4=newcenters;
viscircles(newcenters, newradii,'Color','b');
title(['Frame ',num2str(time),', Slice ', num2str(s)])
newradii=newradii.*(w_cm./ns);
Storage(time).circleradii.slice4=newradii;
CircleArea=pi.*(newradii).^2;
CircleVolume = CircleArea.*longaxisdist
Storage(time).circlearea.slice4=CircleArea;
Storage(time).circlevol.slice4=CircleVolume;
Storage(time).zcenter.slice4=(slice4+slice3)/2; 
%% Area of a Slice
Area = 0;
nonzeros = find(DegreeRadius ~= 0);
count = length(nonzeros);
for i=1:count
    Area = Area + (1/count).*pi.*(DegreeRadius(i).^2);
end
 
Storage(time).newarea.slice4=Area;
Volume = Area.*longaxisdist./count;
MassSlice=density.*Volume;
Storage(time).newvolume.slice4=Volume;
Storage(time).mass.slice4=MassSlice;
end
 
%% Slice5
if s==5
for i=0:thickness
    NewSliceXY = NewHeart(:,:,slice(s) + i);
    ThiccsliceXY = ThiccsliceXY + NewSliceXY;
end
ThiccsliceXY = ThiccsliceXY./ thickness;
 
 
% normalizing averaged slices
ThiccsliceXY=ThiccsliceXY./max(max(ThiccsliceXY));
Filtered2Dxy = zeros(ns);
Filtered2Dxz = zeros(ns);
 
for x = 1:ns
    for y = 1:ns
        %filter xy
        if ThiccsliceXY(x, y) <= 0.36
            Filtered2Dxy(x, y) = 0;
        else
            Filtered2Dxy(x, y) = ThiccsliceXY(x, y);
        end
    end
end
 
figure(7)
imshow(ThiccsliceXY)
title('Thick Slice Before Filtering');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])
 
%% Displaying Thresholded Slices
figure(3)
imshow(Filtered2Dxy);
title('Slice Filtered by Thresholding (Removes Minimum Brightness)');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])
 
%% Filtering 
figure(4);
Ixy=imbinarize(Filtered2Dxy,'global');
edgexy=edge(Ixy,'log');
imshow(Ixy);
title('Filtering Threshold Slice with Otsu Approximation');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])

figure(5)
imshow(edgexy)
title('Edge Detection from Otsu Approximation');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])
 
%% Crop Images to Remove Noise
edgexy = Crop(edgexy, 80, 160, 65, 140, ns); %top bottom left right
 
%% Find Circle
figure(6);
imshow(edgexy)
title('Edge Detection to Generate Circles - Image Has Been Cropped');
[centers,radii] = imfindcircles(edgexy,[8 250]) ;
indradii=find(max(radii)==radii);
for i=1:length(indradii)
    newradii(i)=radii(indradii(i));
    for j=1:2
        newcenters(i,j)= centers(indradii(i),j);
    end
end
%% Area of Circle with 8 Pie Slices (Every 45 Degrees)
DegreeRadius = eightradii(zeros(1, 8), newcenters, edgexy, ns, newradii);
Storage(time).newradii.slice5=DegreeRadius;
%% Starting Storage
Storage(time).circlecenter.slice5=newcenters;
viscircles(newcenters, newradii,'Color','b');
title(['Frame ',num2str(time),', Slice', num2str(s)])
newradii=newradii.*(w_cm./ns);
Storage(time).circleradii.slice5=newradii;
CircleArea=pi.*(newradii).^2;
CircleVolume=CircleArea.*longaxisdist
Storage(time).circlearea.slice5=CircleArea;
Storage(time).circlevol.slice5=CircleVolume;
Storage(time).zcenter.slice5=(slice4+slice5)/2;  
 
%% Area of a Slice
Area = 0;
nonzeros = find(DegreeRadius ~= 0);
count = length(nonzeros);
for i=1:count
    Area = Area + (1/count).*pi.*(DegreeRadius(i).^2);
end
 
Storage(time).newarea.slice5=Area;
Volume = Area.*longaxisdist./count;
MassSlice=density.*Volume;
Storage(time).newvolume.slice5=Volume;
Storage(time).mass.slice5=MassSlice; 
end
 
%% Slice 6
if s==6
for i=0:thickness
    NewSliceXY = NewHeart(:,:,slice(s) + i);
    ThiccsliceXY = ThiccsliceXY + NewSliceXY;
end
ThiccsliceXY = ThiccsliceXY./ thickness;
 
 
% normalizing averaged slices
ThiccsliceXY=ThiccsliceXY./max(max(ThiccsliceXY));
 
%% Attempting to Threshold Filter in 2D
Filtered2Dxy = zeros(ns);
Filtered2Dxz = zeros(ns);
 
for x = 1:ns
    for y = 1:ns
        %filter xy
        if ThiccsliceXY(x, y) <= 0.26
            Filtered2Dxy(x, y) = 0;
        else
            Filtered2Dxy(x, y) = ThiccsliceXY(x, y);
        end
    end
end
 
%% Displaying Thresholded Slices
figure(3)
imshow(Filtered2Dxy);
title('Slice Filtered by Thresholding (Removes Minimum Brightness)');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])
 
%% Filtering 
figure(4);
Ixy=imbinarize(Filtered2Dxy,'global');
edgexy=edge(Ixy,'log');
imshow(Ixy);
title('Filtering Threshold Slice with Otsu Approximation');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])

figure(5)
imshow(edgexy)
title('Edge Detection from Otsu Approximation');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])
 
%% Crop Images to Remove Noise
edgexy = Crop(edgexy, 83, 155, 55, 160, ns);  
 
%% Find Circle
figure(6);
imshow(edgexy)
title('Edge Detection to Generate Circles - Image Has Been Cropped');
[centers,radii] = imfindcircles(edgexy,[8 250]) ;
indradii=find(max(radii)==radii);
for i=1:length(indradii)
    newradii(i)=radii(indradii(i));
    for j=1:2
        newcenters(i,j)= centers(indradii(i),j);
    end
end
 
%% Area of Circle with 8 Pie Slices (Every 45 Degrees)
DegreeRadius = eightradii(zeros(1, 8), newcenters, edgexy, ns, newradii);
Storage(time).newradii.slice6=DegreeRadius;
 
%% Starting Storage
Storage(time).circlecenter.slice6=newcenters;
viscircles(newcenters, newradii,'Color','b');
title(['Frame ',num2str(time),', Slice ', num2str(s)])
newradii=newradii.*(w_cm./ns);
Storage(time).circleradii.slice6=newradii;
CircleArea=pi.*(newradii).^2;
CircleVolume = CircleArea.*longaxisdist
Storage(time).circlearea.slice6=CircleArea;
Storage(time).circlevol.slice6=CircleVolume;
Storage(time).zcenter.slice6=(slice6+slice5)/2;   
%% Area of a Slice
Area = 0;
nonzeros = find(DegreeRadius ~= 0);
count = length(nonzeros);
for i=1:count
    Area = Area + (1/count).*pi.*(DegreeRadius(i).^2);
end
 
Storage(time).newarea.slice6=Area;
Volume = Area.*longaxisdist./count;
MassSlice=density.*Volume;
Storage(time).newvolume.slice6=Volume;
Storage(time).mass.slice6=MassSlice;
end 
 
%% Slice 7
if s==7
for i=0:thickness
    NewSliceXY = NewHeart(:,:,slice(s) + i);
    ThiccsliceXY = ThiccsliceXY + NewSliceXY;
end
ThiccsliceXY = ThiccsliceXY./ thickness;
 
 
% normalizing averaged slices
ThiccsliceXY=ThiccsliceXY./max(max(ThiccsliceXY));
Filtered2Dxy = zeros(ns);
Filtered2Dxz = zeros(ns);
 
for x = 1:ns
    for y = 1:ns
        %filter xy
        if ThiccsliceXY(x, y) <= 0.39
            Filtered2Dxy(x, y) = 0;
        else
            Filtered2Dxy(x, y) = ThiccsliceXY(x, y);
        end
    end
end
 
figure(7)
imshow(ThiccsliceXY)
title('Thick Slice Before Filtering');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])
 
%% Displaying Thresholded Slices
figure(3)
imshow(Filtered2Dxy);
title('Slice Filtered by Thresholding (Removes Minimum Brightness)');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])
 
%% Filtering 
figure(4);
Ixy=imbinarize(Filtered2Dxy,'global');
edgexy=edge(Ixy,'log');
imshow(Ixy);
title('Filtering Threshold Slice with Otsu Approximation');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])

figure(5)
imshow(edgexy)
title('Edge Detection from Otsu Approximation');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])
 
%% Crop Images to Remove Noise
edgexy = Crop(edgexy, 90, 150, 60, 140, ns);
 
%% Find Circle
figure(6);
imshow(edgexy)
title('Edge Detection to Generate Circles - Image Has Been Cropped');
[centers,radii] = imfindcircles(edgexy,[8 250]) ;
indradii=find(max(radii)==radii);
for i=1:length(indradii)
    newradii(i)=radii(indradii(i));
    for j=1:2
        newcenters(i,j)= centers(indradii(i),j);
    end
end
 
%% Area of Circle with 8 Pie Slices (Every 45 Degrees)
DegreeRadius = eightradii(zeros(1, 8), newcenters, edgexy, ns, newradii);
Storage(time).newradii.slice7=DegreeRadius;
 
%% Starting Storage
Storage(time).circlecenter.slice7=newcenters;
viscircles(newcenters, newradii,'Color','b');
title(['Frame ',num2str(time),', Slice ', num2str(s)])
newradii=newradii.*(w_cm./ns);
Storage(time).circleradii.slice7=newradii;
CircleArea=pi.*(newradii).^2;
CircleVolume=CircleArea.*longaxisdist
Storage(time).circlearea.slice7=CircleArea;
Storage(time).circlevol.slice7=CircleVolume;
Storage(time).zcenter.slice7=(slice6+slice7)/2;    
%% Area of a Slice
Area = 0;
nonzeros = find(DegreeRadius ~= 0);
count = length(nonzeros);
for i=1:count
    Area = Area + (1/count).*pi.*(DegreeRadius(i).^2);
end
 
Storage(time).newarea.slice7=Area;
Volume = Area.*longaxisdist./count;
MassSlice=density.*Volume;
Storage(time).newvolume.slice7=Volume;
Storage(time).mass.slice7=MassSlice;
end
 
%% Slice 8 
if s==8
for i=0:thickness
    NewSliceXY = NewHeart(:,:,slice(s) + i);
    ThiccsliceXY = ThiccsliceXY + NewSliceXY;
end
ThiccsliceXY = ThiccsliceXY./ thickness;
 
 
% normalizing averaged slices
ThiccsliceXY=ThiccsliceXY./max(max(ThiccsliceXY));
Filtered2Dxy = zeros(ns);
Filtered2Dxz = zeros(ns);
 
for x = 1:ns
    for y = 1:ns
        %filter xy
        if ThiccsliceXY(x, y) <= 0.45
            Filtered2Dxy(x, y) = 0;
        else
            Filtered2Dxy(x, y) = ThiccsliceXY(x, y);
        end
    end
end
 
figure(7)
imshow(ThiccsliceXY)
title('Thick Slice Before Filtering');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])
 
%% Displaying Thresholded Slices
figure(3)
imshow(Filtered2Dxy);
title('Slice Filtered by Thresholding (Removes Minimum Brightness)');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])
 
%% Filtering 
figure(4);
Ixy=imbinarize(Filtered2Dxy,'global');
edgexy=edge(Ixy,'log');
imshow(Ixy);
title('Filtering Threshold Slice with Otsu Approximation');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])

figure(5)
imshow(edgexy)
title('Edge Detection from Otsu Approximation');
xlabel(['Frame ',num2str(time),', Slice ', num2str(s)])
 
 
%% Crop Images to Remove Noise
edgexy = Crop(edgexy, 96 , 140, 74, 140, ns);
 
%% Find Circle
figure(6);
imshow(edgexy)
title('Edge Detection to Generate Circles - Image Has Been Cropped');
[centers,radii] = imfindcircles(edgexy,[8 250]) ;
indradii=find(max(radii)==radii);
for i=1:length(indradii)
    newradii(i)=radii(indradii(i));
    for j=1:2
        newcenters(i,j)= centers(indradii(i),j);
    end
end
 
%% Area of Circle with 8 Pie Slices (Every 45 Degrees)
DegreeRadius = eightradii(zeros(1, 8), newcenters, edgexy, ns, newradii);
Storage(time).newradii.slice8=DegreeRadius;
 
%% Starting Storage
Storage(time).circlecenter.slice8=newcenters;
viscircles(newcenters, newradii,'Color','b');
title(['Frame ',num2str(time),', Slice ', num2str(s)])
newradii=newradii.*(w_cm./ns);
Storage(time).circleradii.slice8=newradii;
CircleArea=pi.*(newradii).^2;
CircleVolume=CircleArea.*longaxisdist
Storage(time).circlearea.slice8=CircleArea;
Storage(time).circlevol.slice8=CircleVolume;
Storage(time).zcenter.slice8=(slice8+slice7)/2; 
%% Area of a Slice
Area = 0;
nonzeros = find(DegreeRadius ~= 0);
count = length(nonzeros);
for i=1:count
    Area = Area + (1/count).*pi.*(DegreeRadius(i).^2);
end
 
Storage(time).newarea.slice8=Area;
Volume = Area.*longaxisdist./count;
MassSlice=density.*Volume;
Storage(time).newvolume.slice8=Volume;
Storage(time).mass.slice8=MassSlice;
end 
 
%% Pause to view figures before going to next time frame
pause
close all
clear centers radii newradii newcenters indradii edgexy X Y X0 Y0 Point1 Point2 Point3 Point4 Point5 Point6 Point7 Point8 DegreeRadius
end
%% Extract Data
% Center Coordinates
% x and y coordinate averaging
Storage(time).avgxy=1/8.*(Storage(time).circlecenter.slice1 + ...
    Storage(time).circlecenter.slice2 + Storage(time).circlecenter.slice3+ ...
    Storage(time).circlecenter.slice4 + Storage(time).circlecenter.slice5+ ...
    Storage(time).circlecenter.slice6 + Storage(time).circlecenter.slice7+ ...
    Storage(time).circlecenter.slice8);
 
% Circle Volume
Storage(time).avgcirclevolume=1/8.*(Storage(time).circlevol.slice1 + ...
    Storage(time).circlevol.slice2 + Storage(time).circlevol.slice3 + ...
    Storage(time).circlevol.slice4 + Storage(time).circlevol.slice5 + ...
    Storage(time).circlevol.slice6 + Storage(time).circlevol.slice7 + ...
    Storage(time).circlevol.slice8);
 
% New Volume
Storage(time).avgvolume=(Storage(time).newvolume.slice1 + ...
    Storage(time).newvolume.slice2 + Storage(time).newvolume.slice3 + ...
    Storage(time).newvolume.slice4 + Storage(time).newvolume.slice5 + ...
    Storage(time).newvolume.slice6 + Storage(time).newvolume.slice7 + ...
    Storage(time).newvolume.slice8);
% Total Mass
Storage(time).totalmass=Storage(time).avgvolume.*density;

% Center of Mass Equation for x
Storage(time).xcom=(1/Storage(time).totalmass).*(Storage(time).circlecenter.slice1(1).*Storage(time).mass.slice1 + ...
    Storage(time).circlecenter.slice2(1).*Storage(time).mass.slice2 + ...
    Storage(time).circlecenter.slice3(1).*Storage(time).mass.slice3 + ...
    Storage(time).circlecenter.slice4(1).*Storage(time).mass.slice4 + ...
    Storage(time).circlecenter.slice5(1).*Storage(time).mass.slice5 + ...
    Storage(time).circlecenter.slice6(1).*Storage(time).mass.slice6 + ...
    Storage(time).circlecenter.slice7(1).*Storage(time).mass.slice7 + ...
    Storage(time).circlecenter.slice8(1).*Storage(time).mass.slice8);
% Center of Mass Equation for y
Storage(time).ycom=(1/Storage(time).totalmass).*(Storage(time).circlecenter.slice1(2).*Storage(time).mass.slice1 + ...
    Storage(time).circlecenter.slice2(2).*Storage(time).mass.slice2 + ...
    Storage(time).circlecenter.slice3(2).*Storage(time).mass.slice3 + ...
    Storage(time).circlecenter.slice4(2).*Storage(time).mass.slice4 + ...
    Storage(time).circlecenter.slice5(2).*Storage(time).mass.slice5 + ...
    Storage(time).circlecenter.slice6(2).*Storage(time).mass.slice6 + ...
    Storage(time).circlecenter.slice7(2).*Storage(time).mass.slice7 + ...
    Storage(time).circlecenter.slice8(2).*Storage(time).mass.slice8);
% Center of Mass Equation for z
Storage(time).zcom=(1/Storage(time).totalmass).*(Storage(time).zcenter.slice1.*Storage(time).mass.slice1 + ...
    Storage(time).zcenter.slice2.*Storage(time).mass.slice2 + ...
    Storage(time).zcenter.slice3.*Storage(time).mass.slice3 + ...
    Storage(time).zcenter.slice4.*Storage(time).mass.slice4 + ...
    Storage(time).zcenter.slice5.*Storage(time).mass.slice5 + ...
    Storage(time).zcenter.slice6.*Storage(time).mass.slice6 + ...
    Storage(time).zcenter.slice7.*Storage(time).mass.slice7 + ...
    Storage(time).zcenter.slice8.*Storage(time).mass.slice8);

end
