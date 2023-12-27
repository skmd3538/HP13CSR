%returns intermediate steps
function [subsampled, interpolatedBicubic, lowResImage] = prepareLowRes(input, tgtSiz)

%simple lowpass/blur filter, all 1s/(size*size), keep size at 3
filterSize = 3;
blurFilter = ones(filterSize,filterSize) * 1/(filterSize*filterSize);
blurredImage = im2double(imfilter(input, blurFilter));

%subsample and interpolate
inputSiz = size(input);
if ~exist('tgtSiz','var')
    tgtSiz = inputSiz / 2; 
end


% subsampled = blurredImage(1:2:end,1:2:end,:);
skip = floor(inputSiz(1)/tgtSiz(1));
idx = 1:skip:size(blurredImage, 1);
subsampled = blurredImage(idx,idx,:);
% subsampled = imresize(blurredImage, tgtSiz, 'nearest');
interpolatedBicubic = imresize(subsampled, inputSiz(1:2), 'bicubic');

%simple high pass filter, all pass - low pass
highPassFilter = [0 0 0; 0 1 0; 0 0 0] - blurFilter;
filteredImage = imfilter(interpolatedBicubic, highPassFilter);

lowResImage = filteredImage;
end




