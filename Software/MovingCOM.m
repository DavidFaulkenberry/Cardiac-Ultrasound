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
 

%% Moving COM Plot
%AllPoints = zeros(1, 25);
AllPoints(time) = time;
Points = round(AllPoints);


figure(1)
X = [98.45553347
98.98399125
98.57407344
98.17912674
98.69303169
98.78967727
98.67718578
98.84229919
98.94955871
99.36191746
99.76104862
99.44183138
98.10580693
99.95975913
95.7975048
99.93081445
99.57551483
98.92853444
100.2199258
98.98122908
99.10826578
99.07019694
99.00491866
99.79076719
100.3636363];% Average X COM for each Frame

Y = [114.2094249
114.8783686
114.6148522
114.8760814
114.2362137
115.2627392
114.7403606
115.0244461
115.3411608
116.6007166
115.8968634
117.4513095
116.4843968
117.9864494
114.6874184
115.3922918
116.5657695
116.8805006
114.8711252
116.3286877
116.2816841
116.2544812
116.2929633
115.88823
115.8231863]; % Average Y COM for each Frame

Z = [91.92903982
93.21550719
94.86660431
93.7687467
91.5175285
92.63403888
91.54884928
92.90826604
92.58544409
93.34433151
91.80715513
93.3566173
95.35980202
93.59753474
94.38630762
96.49160434
91.22665741
90.7601913
87.74076891
93.29156455
92.80709703
92.80182033
92.90181452
90.30289184
93.81172691]; % Average Z COM for each Frame
   
    subplot(1, 2, 1)
    imshow(xyview);
    hold on
    plot(X(time), Y(time), 'r.');
    xlabel('XY View');
    
    subplot(1, 2, 2)
    imshow(xzview);
    hold on;
    plot(X(time), Z(time), 'r.');
    xlabel('XZ View');
    suptitle('Center of Mass Movement Over Time');
    pause(0.01)

    % Input time to get singular point on each frame, input Points to get
    % compounding image of points over time
end