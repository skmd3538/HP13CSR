function [keys, values, imgSiz, tgtSiz, params] = train(inputImagePath, params)

%read in input image
[~, ~, fExt] = fileparts(inputImagePath);
if(isdicom(inputImagePath))
    inputImage = dicomread(inputImagePath);
elseif strcmp(fExt, '.mat')
    inputImage = load(inputImagePath);
    inputImage = abs(struct2array(inputImage(1)));
else
    inputImage = imread(inputImagePath);
end

inputImage = mat2gray(im2double(inputImage));


if(isfield(params, 'tgtSiz'))
    tgtSiz = params.tgtSiz;
elseif(isfield(params, 'scale'))
    params.tgtSiz = size(inputImage) * params.scale;
    tgtSiz = params.tgtSiz;
else
    params.tgtSiz = size(inputImage);
    tgtSiz = params.tgtSiz;
end


if(size(inputImage, 3) < 3)
    inputImage = repmat(inputImage, [1 1 3]);
end


%perform preprocessing steps, blur/downsample/interpolate/highpass
if(isfield(params, 'externalFile'))
    [subSampled, interpolatedSubsampled, lowResImage] = prepareLowResExternal(params);
else
    [subSampled, interpolatedSubsampled, lowResImage] = prepareLowRes(inputImage, tgtSiz);
end

highResImage = im2double(inputImage) - im2double(interpolatedSubsampled);

if(~isfield(params, 'LRPatchSiz'))
    params.LRPatchSiz = 7;
end

if(~isfield(params, 'HRPatchSiz'))
    params.HRPatchSiz = 5;
end

%make lowResPatches
[lowResPatches, lowResRows, lowResCols] = make2dPatchMap(lowResImage, params.LRPatchSiz);
%DO NOT USE HIGH RES DIMENSIONS
[highResPatches , numRows, numCols] = make2dPatchMap(highResImage, params.HRPatchSiz);
   
%create key/values based on patches
alpha = params.alpha;
[keys, values] = createVectors(alpha, lowResPatches, highResPatches, lowResRows, lowResCols);

imgSiz = size(subSampled);

end