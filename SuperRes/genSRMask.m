% Mask generation
%   Mask = genSRMask(iMag, voxel_size)
% 
%   output
%   Mask - the biggest contiguous region that has decent SNR
%
%   input
%   iMag - the magnitude MR image
%   voxel_size - the size of a voxel
%


function Mask = genSRMask(iMag)
% if nargin < 3; erosion = 0.999; end
    %matrix_size = size(iMag);
    m1 = iMag>(0.05*max(iMag(:)));           % simple threshold
%     m1 = SMV(m,matrix_size, voxel_size, 10)>erosion;   % erode the boundary by 10mm
%     l = bwlabeln(m1,6);                       % find the biggest contiguous region
%     Mask = (l==mode(l(l~=0)));
%     Mask1 = SMV(Mask, matrix_size, voxel_size, 10)>0.001; % restore the enrosion by 4mm
    Mask = m1;
end
