function [subsampled, interpolated, superResImage] = rebuildImage(kdTree, subSampledInput, values, params)

%recreate the low res image that you will be adding the difference to
tgtSiz = params.tgtSiz;
backImage = im2double(imresize(subSampledInput, tgtSiz(1:2), params.interpMethod));
% backImage = im2double(imresize(subSampledInput, 2.0, 'bicubic'));
interpolatedSubSampledInput = backImage;
[r, c, d] = size(backImage);

%obtain high frequencies from the low res image
filterSize = 3;
blurFilter = ones(filterSize,filterSize) * 1/(filterSize*filterSize);

highPassFilter = [0 0 0; 0 1 0; 0 0 0] - blurFilter;
filteredImage = imfilter(backImage, highPassFilter);


%diff difference between low res and high passed version of the low res image
%it is added back to reconstructed image at the end
diff = backImage - filteredImage;
backImage = filteredImage;


% [patches, numRows, numCols] = make2dPatchMap(backImage, params.LRPatchSiz);


%initialize reconstructed image
reconstructedImage = zeros(r,c,d);

alpha = params.alpha;
lrLen = params.LRPatchSiz - 1;
hrLen = params.HRPatchSiz;
patchLen = (params.LRPatchSiz  * params.LRPatchSiz  + ...
    (params.HRPatchSiz + params.HRPatchSiz - 1));
endR = r-lrLen;
endC = c-lrLen;

%while (i + 6 < r) 
for i = 1:hrLen-1:endR
    %print current percentage done
%     percentage = i/endR*100
    for j = 1:hrLen-1:endC

    %get corresponding low res patch
    patch49 = backImage(i:i+lrLen, j:j+lrLen, :);
    
    %obtain current high res patch being built for extracting
    %9 upper left values
    patch25 = reconstructedImage((i+1):(i+hrLen), (j+1):(j+hrLen), :);
      
     
    %concatenate patch so has 58 elements and becomes key
    patch58 = zeros(patchLen,3);
    for k = 1:3
        patch9 = [];
        patch9 = [patch25(1,:,k) patch25(2:end,1,k)'];
        patch9 = patch9*alpha;
        patch58(:,k) = [reshape(patch49(:,:,k), 1, []) patch9];
    end        
    patch58 = reshape(patch58, 1, 3*patchLen);
    
    %undo contrast normalize
    scale = getContrastNormalizeScale(patch58);
    patch58 = patch58 / scale;

    %obtain the value
    index = knnsearch(kdTree, patch58);
    value(:,:,:) = values(index,:,:,:);

    %add high res patch to the image being reconstructed
    reconstructedImage((i+1):(i+hrLen), (j+1):(j+hrLen), :) = value*scale;
      
    end
end

output = reconstructedImage+backImage+diff;

superResImage = output;


subsampled = subSampledInput;
interpolated = interpolatedSubSampledInput;


end