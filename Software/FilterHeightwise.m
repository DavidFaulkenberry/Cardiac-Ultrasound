function vol_filtered = FilterHeightwise(vol)

% Blackman window size N is the filter
 N = 21;
 % 1D filter as a 3D point spread function
 filter = zeros(1,N,1);
 filter(1,:,1) = blackman(N);

 % PSF convolution in 3D with the 1D filter (nothing will happen in the
 % [width,depth] plane as the convolution is along every height vector as
 % if its own 1D signal), and then cut off the extra N/2 added on each
 % side of the height dimension to maintain unchanged height dimension
 % size (this should have a negligible effect on our analysis).
 vol_filtered = convn(vol,filter,'same');
end