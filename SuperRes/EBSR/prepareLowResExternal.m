%returns intermediate steps
function [subsampled, interpolatedBicubic, lowResImage] = prepareLowResExternal(params)

inputImage = params.externalFileName;
%read in input image
[~, ~, fExt] = fileparts(inputImage);
if(isdicom(inputImage))
    input = dicomread(inputImage);
elseif strcmp(fExt, '.mat')
    input = load(inputImage);
    input = abs(struct2array(input(1)));
else
    input = imread(inputImage);
end

% if(params.grayscale ~= 1)
   input = mat2gray(im2double(input));
% end



%simple lowpass/blur filter, all 1s/(size*size), keep size at 3
filterSize = 3;
blurFilter = ones(filterSize,filterSize) * 1/(filterSize*filterSize);

%subsample and interpolate
% inputSiz = size(input);


subsampled = input;
% subsampled = imresize(blurredImage, tgtSiz, 'nearest');
% interpolatedBicubic = imresize(subsampled, inputSiz(1:2)*params.scale, params.interpMethod);
interpolatedBicubic = imresize(subsampled, params.tgtSiz, params.interpMethod);

%simple high pass filter, all pass - low pass
% if(params.scale ~= 1)
    highPassFilter = [0 0 0; 0 1 0; 0 0 0] - blurFilter;
    filteredImage = imfilter(interpolatedBicubic, highPassFilter);
% else
%     filteredImage = interpolatedBicubic;%just denoising when no resizing
% end

lowResImage = filteredImage;
end




