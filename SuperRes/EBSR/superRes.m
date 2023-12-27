function [subsampled, interpolated, MainSR] = superRes(fileName, inputImageName, params)
% 131-139 dicomwrite
    
load(fileName, 'keys', 'values', 'LRPatchSiz', 'HRPatchSiz');

bucketSize = params.bucketSize;

kdTree = KDTreeSearcher(keys, 'BucketSize', bucketSize);
sz3 = 1;

if(~isnumeric(inputImageName))
    [~, ~, fExt] = fileparts(inputImageName);
end
%convert original image to double version
if(isnumeric(inputImageName))
    subsampledInput = inputImageName;
     params.sz = size(subsampledInput);
     subsampledInput = reshape(subsampledInput, ...
         params.sz(1), params.sz(2), []);
     sz3 = prod(params.sz(3:end));
elseif(isdicom(inputImageName))
    subsampledInput = dicomread(inputImageName);
     params.sz = size(subsampledInput);
elseif strcmp(fExt, '.mat')
     subsampledInput = load(inputImageName);  
     subsampledInput = abs(struct2array(subsampledInput(1)));

     nd = ndims(subsampledInput);
     params.sz = size(subsampledInput);
     if nd > 2
         if(isfield(params, 'tp'))
             subsampledInput = reshape(subsampledInput, [], params.sz(end));
             subsampledInput = subsampledInput(:, 1:params.tp);
         end
         subsampledInput = reshape(subsampledInput, params.sz(1), params.sz(2), []);
         sz3 = size(subsampledInput,3);
     end
else
    subsampledInput = imread(inputImageName);
    params.sz = size(subsampledInput);
end

%subsampledInput = mat2gray(im2double(subsampledInput));


inputSiz = size(subsampledInput);
if(isfield(params, 'tgtSiz'))
    tgtSiz = params.tgtSiz;
else
    tgtSiz = inputSiz(1:2) * 2; 
end

params.tgtSiz = tgtSiz;
params.inputSiz = inputSiz(1:2);

params.MainSR = zeros(tgtSiz(1), tgtSiz(2), prod(params.sz(3:end)));
params.interP = zeros(tgtSiz(1), tgtSiz(2), prod(params.sz(3:end)));
params.subsampledInput = subsampledInput;

% cleanupObj = onCleanup(@()cleanUpFun(params));
  
params.HRPatchSiz = HRPatchSiz;
params.LRPatchSiz = LRPatchSiz;

% cleanupObj = onCleanup(@()cleanUpFun(params));


for i = 1:sz3
  
    fprintf('Processing frame %d of %d\n', i, sz3);
%      input = subsampledInput(:,:,i);
    %Scale all images to range [0 1] to be compatible with weights
    input = mat2gray(subsampledInput(:,:,i));
    if(size(input, 3) < 3)
        input = repmat(input, [1 1 3]);
    end
    
    
    %low frequency and subsampled version of original image to use
    %for the actually super resolution and to compare results

    
    [subsampled, interpolated, superResImage] = rebuildImage(kdTree, input, values, params);
    
    %rescale image grayscale values to original range 
    if(params.rescale)
        res = rescale(superResImage(:,:,1), min(min(subsampledInput(:,:,i))) ...
            ,  max(max(subsampledInput(:,:,i))));
    else
        res = superResImage(:,:,1);
    end
    params.MainSR(:,:,i) = res;
    
    res = rescale(interpolated(:,:,1), min(min(subsampledInput(:,:,i))) ...
         ,  max(max(subsampledInput(:,:,i))));
  %  res = interpolated(:,:,1);
    params.interP(:,:,i) = res;
    
%     if(isdicom(inputImageName))

   
end


saveData(params);

MainSR = params.MainSR;


end



% fires when main function terminates
function saveData(param)
savePath = sprintf('%s_EBSR_%dto%d', param.savePath, param.inputSiz(1), param.tgtSiz(1));

if(~exist(savePath, 'dir'))
    mkdir(savePath);
end

sz = param.sz;
if(isfield(param, 'tp'))
    sz(end) = param.tp;
end
sz(1:2) = param.tgtSiz(1:2);
EBSR = reshape(param.MainSR, sz);
interP = reshape(param.interP, sz);
NI = imresize(param.subsampledInput, sz(1:2), 'nearest');
NI = double(reshape(NI, sz));
BCI = imresize(param.subsampledInput, sz(1:2), 'bicubic');
BCI = double(reshape(BCI, sz));
SI = imresize(param.subsampledInput, sz(1:2), 'lanczos2');
SI = double(reshape(SI, sz));
% SI = sincinterp (param.subsampledInput, sz(1),sz(2)); %sinc interpolation
% SI = reshape(SI, sz);
save(fullfile(savePath, 'superRes.mat'), 'EBSR', 'SI', 'NI', 'BCI', 'interP');
[p,n,s]=fileparts(savePath);
save(fullfile(savePath, sprintf('%s.mat',n)), 'EBSR')



if(isfield(param, 'DicomTemplate'))
    
    dicomTemplate = param.DicomTemplate;
  

    dicFiles = dir(fullfile(dicomTemplate, '*.dcm'));
    dic_info = dicominfo(fullfile(dicomTemplate,dicFiles(1).name));
    
    dic_info.PixelSpacing(1) = dic_info.PixelSpacing(1) * double(dic_info.Width)/param.tgtSiz(1);
    dic_info.PixelSpacing(2) = dic_info.PixelSpacing(2) * double(dic_info.Height)/param.tgtSiz(2);
    dic_info.Width = param.tgtSiz(1);
    dic_info.Height = param.tgtSiz(2);
    dic_info.Rows = param.tgtSiz(1);
    dic_info.Columns = param.tgtSiz(2);
    if(isfield(param, 'msk_header') && contains(param.msk_header.sequence, 'CSI'))
%         dic_info.ImagePositionPatient = param.msk_header.toplc;
    end
  
    metFolder = param.LRName;
    for i = 1:length(param.LRName)
        metFolder{i} = fullfile(savePath, param.LRName{i});
        if( 7 ~= exist(metFolder{i}, 'dir'))
            mkdir(metFolder{i});
        end
    end
    
    
    met_map = EBSR;
    sz = size(met_map);
    %If only 2D, set slice, time and num mets dims to 1
    if(length(sz) < 5)
        for i = length(sz)+1:5
            sz(i) = 1;
        end
    end
    
    if(isfield(dic_info, 'ImagesInAcquisition'))
        nImgs = dic_info.ImagesInAcquisition;
    else
        nImgs = 1;
    end
    imgSave = reshape(met_map, sz(1), sz(2), [], sz(5));
    nImgsSave = size(imgSave, 3);
    oneFrame = nImgs/nImgsSave;
    
    
    
    for i = 1:size(imgSave, 4) %metabolites
        for j = 1:size(imgSave, 3) %slices x timepoints
            img = imgSave(:,:,j, i);
            dic_info2 = dicominfo(fullfile(dicomTemplate,dicFiles(oneFrame*j).name));
            if(isfield(dic_info2, 'SmallestImagePixelValue') && isfield(dic_info2, 'LargestImagePixelValue'))
                img = rescale(img, dic_info.SmallestImagePixelValue, ...
                    dic_info2.LargestImagePixelValue );
            end
            dic_info.InstanceNumber = j;
            dic_info.AcquisitionNumber = j;
            dic_info.ImagesInAcquisition = nImgsSave;
            dic_info.SeriesDescription = sprintf('%s_%d', param.LRName{i}, j);
            dicomwrite(uint16(img), fullfile(metFolder{i}, sprintf('%s_%04d_%dx%d.dcm', ...
                param.LRName{i}, j, param.tgtSiz(1), param.tgtSiz(2))), dic_info, 'CreateMode', 'copy');
        end
    end
    
    
    
end



if isfield(param, 'TruthImage')
    
    [~, ~, fEx] = fileparts(param.TruthImage);
    if(isdicom(param.TruthImage))
        TI  = im2double(dicomread(param.TruthImage));
    elseif strcmp(fEx, '.mat')
        TI  = load(param.TruthImage);
        TI  = im2double(abs(struct2array(TI (1))));
    else
        TI = im2double(imread(param.TruthImage));
    end
    
     TI = mat2gray(TI);

    cmap = colormap(jet(64));
    fontSize = 10;    
    
    sz3 = size(EBSR, 3);
    
    for sl = 1:sz3
        
    ti = double(TI(:,:,sl));
    si = mat2gray(double(SI(:,:,sl)));
    ni = mat2gray(double(NI(:,:,sl)));
    ebsr = double(EBSR(:,:,sl));
    bci = mat2gray(double(BCI(:,:,sl)));
    
    subplot(sz3,5,sub2ind([5,sz3],1,sl)),imshow(mat2gray(ti), 'DisplayRange', [],'ColorMap', cmap),title('Truth','FontSize', fontSize);hold  on
    subplot(sz3,5,sub2ind([5,sz3],2,sl)),imshow(mat2gray(si), 'DisplayRange', [],'ColorMap', cmap),title('Lanczos2','FontSize', fontSize);
    mseSI = immse(ti,si);
    ssimSI = ssim(ti,si);
    psnrSI = psnr(ti, si);
    xlabel(sprintf('PSNR:%.03f, SSIM:%.03f', psnrSI, ssimSI), 'FontSize', fontSize);
    subplot(sz3,5,sub2ind([5,sz3],3,sl)),imshow(mat2gray(ni), 'DisplayRange', [],'ColorMap', cmap),title('Nearest','FontSize', fontSize);
    mseNI = immse(ti,ni);
    ssimNI = ssim(ti,ni);
    psnrNI = psnr(ti, ni);
    xlabel(sprintf('PSNR:%.01f, SSIM:%.03f', psnrNI, ssimNI), 'FontSize', fontSize);
    subplot(sz3,5,sub2ind([5,sz3],4,sl)),imshow(mat2gray(bci), 'DisplayRange', [],'ColorMap', cmap),title('Bicubic','FontSize', fontSize);
    mseBCI = immse(ti,bci);
    ssimBCI = ssim(ti,bci);
    psnrBCI = psnr(ti, bci);
    xlabel(sprintf('PSNR:%.01f, SSIM:%.03f', psnrBCI, ssimBCI), 'FontSize', fontSize);
    subplot(sz3,5,sub2ind([5,sz3],5,sl)),imshow(mat2gray(ebsr), 'DisplayRange', [],'ColorMap', cmap),title('IQT','FontSize', fontSize);hold off
    mseSR = immse(ti,ebsr);
    ssimSR = ssim(ti,ebsr);
    psnrSR = psnr(ti, ebsr);
    xlabel(sprintf('PSNR:%.03f, SSIM:%.03f', psnrSR, ssimSR), 'FontSize', fontSize);

    end

    save(fullfile(savePath, 'superRes.mat'), 'TI', 'SI', 'EBSR', 'NI', 'BCI', 'interP');
    save(fullfile(savePath, 'EBSR.mat'), 'EBSR');

end
end


function cleanUpFun(param)

 saveData(param);
end