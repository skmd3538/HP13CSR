function [subsampled, interpolated, MainSR] = superRes(fileName, inputImageName, params)
% 131-139 dicomwrite
    
load(fileName, 'keys', 'values', 'LRPatchSiz', 'HRPatchSiz');

bucketSize = params.bucketSize;

kdTree = KDTreeSearcher(keys, 'BucketSize', bucketSize);
sz3 = 1;

[~, ~, fExt] = fileparts(inputImageName);
%convert original image to double version
if(isdicom(inputImageName))
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

% if(params.grayscale ~= 1)
%     subsampledInput = mat2gray(im2double(subsampledInput));
% end


inputSiz = size(subsampledInput);
if(~isfield(params, 'tgtSiz'))
     params.tgtSiz = inputSiz(1:2) * 2;
end

params.inputSiz = inputSiz(1:2);

params.MainSR = zeros(params.tgtSiz(1), params.tgtSiz(2), prod(params.sz(3:end)));
params.interP = zeros(params.tgtSiz(1), params.tgtSiz(2), prod(params.sz(3:end)));
params.subsampledInput = subsampledInput;

% cleanupObj = onCleanup(@()cleanUpFun(params));
  
params.HRPatchSiz = HRPatchSiz;
params.LRPatchSiz = LRPatchSiz;

% cleanupObj = onCleanup(@()cleanUpFun(params));


for i = 1:sz3
  
    fprintf('Processing frame %d of %d\n', i, sz3);
     input = subsampledInput(:,:,i);
     
     if(params.grayscale ~= 1)
         input = mat2gray(im2double(input));
     end
     
    if(size(input, 3) < 3)
        input = repmat(input, [1 1 3]);
    end
    
    
    %low frequency and subsampled version of original image to use
    %for the actually super resolution and to compare results

    
    [subsampled, interpolated, superResImage] = rebuildImage(kdTree, input, values, params);
    
    params.MainSR(:,:,i) = superResImage(:,:,1);
    params.interP(:,:,i) = interpolated(:,:,1);
    
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
BLI = imresize(param.subsampledInput, sz(1:2), 'bilinear');
BLI = double(reshape(BLI, sz));
SI = imresize(param.subsampledInput, sz(1:2), 'lanczos3');
SI = double(reshape(SI, sz));
% SI = sincinterp (param.subsampledInput, sz(1),sz(2)); %sinc interpolation
% SI = reshape(SI, sz);
save(fullfile(savePath, 'superRes.mat'), 'EBSR', 'SI', 'NI', 'BLI', 'interP');
    save(fullfile(savePath, 'EBSR.mat'), 'EBSR');
fprintf('Results saved to %s\n', savePath);


if(isfield(param, 'DicomTemplate'))
    
    dicomTemplate = param.DicomTemplate;
    [~, f, ~] = fileparts(param.DicomTemplate);
    dic_info = dicominfo(dicomTemplate);
    dic_info.PixelSpacing(1) = dic_info.PixelSpacing(1) * double(dic_info.Width)/param.tgtSiz(1);
    dic_info.PixelSpacing(2) = dic_info.PixelSpacing(2) * double(dic_info.Height)/param.tgtSiz(2);
    dic_info.Width = param.tgtSiz(1);
    dic_info.Height = param.tgtSiz(2);
    dic_info.Rows = param.tgtSiz(1);
    dic_info.Columns = param.tgtSiz(2);
    dic_info.SeriesDescription = sprintf('EBSR%s%d', param.LRName, param.tgtSiz(1));
%    if(~strcmp(f, param.LRName))
     if(false)
%         img = mat2gray(rot90(EBSR,1));
img =  rot90(EBSR,-2);
% img =  EBSR;
        dicomwrite(img, fullfile(savePath, sprintf('%s_SR%d.dcm', param.LRName, param.tgtSiz(1))), dic_info);
    else
%         img =  mat2gray(EBSR);
         img =  EBSR;
                 dicomwrite(img, fullfile(savePath, sprintf('%s_SR%d.dcm', param.LRName, param.tgtSiz(1))), dic_info);
%         dicomwrite(img, fullfile(savePath, sprintf('%s_SR%d.dcm', param.LRName, param.scale)));
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
    
    if(param.grayscale ~= 1)
        TI = mat2gray(TI);
    end

    cmap = colormap(jet(64));
    fontSize = 20;    
    
    sz3 = size(EBSR, 3);
    
    for sl = 1:sz3
        
    ti = TI(:,:,sl);
    si = SI(:,:,sl);
    ni = NI(:,:,sl);
    ebsr = EBSR(:,:,sl);
    intp = interP(:,:,sl);
    
    
    if(param.grayscale ~= 1)
        ti = mat2gray(ti);
        si = mat2gray(si);
        ni = mat2gray(ni);
        ebsr = mat2gray(ebsr);
        intp = mat2gray(intp);
    end
    
    subplot(sz3,5,sub2ind([5,sz3],1,sl)),imshow(ti, 'DisplayRange', [],'ColorMap', cmap),title('Truth','FontSize', fontSize);hold  on
    subplot(sz3,5,sub2ind([5,sz3],2,sl)),imshow(si, 'DisplayRange', [],'ColorMap', cmap),title('Zerofill','FontSize', fontSize);
    mseSI = immse(ti,si);
    ssimSI = ssim(ti,si);
    psnrSI = psnr(ti, si);
    xlabel(sprintf('PSNR:%.03f, SSIM:%.03f', psnrSI, ssimSI), 'FontSize', fontSize);
    subplot(sz3,5,sub2ind([5,sz3],3,sl)),imshow(ni, 'DisplayRange', [],'ColorMap', cmap),title('Nearest','FontSize', fontSize);
    mseNI = immse(ti,ni);
    ssimNI = ssim(ti,ni);
    psnrNI = psnr(ti, ni);
    xlabel(sprintf('PSNR:%.01f, SSIM:%.03f', psnrNI, ssimNI), 'FontSize', fontSize);
    subplot(sz3,5,sub2ind([5,sz3],4,sl)),imshow(intp, 'DisplayRange', [],'ColorMap', cmap),title('Bicubic','FontSize', fontSize);
    mseBCI = immse(ti,intp);
    ssimBCI = ssim(ti,intp);
    psnrBCI = psnr(ti, intp);
    xlabel(sprintf('PSNR:%.01f, SSIM:%.03f', psnrBCI, ssimBCI), 'FontSize', fontSize);
    subplot(sz3,5,sub2ind([5,sz3],5,sl)),imshow(ebsr, 'DisplayRange', [],'ColorMap', cmap),title('EBSR','FontSize', fontSize);hold off
    mseSR = immse(ti,ebsr);
    ssimSR = ssim(ti,ebsr);
    psnrSR = psnr(ti, ebsr);
    xlabel(sprintf('PSNR:%.03f, SSIM:%.03f', psnrSR, ssimSR), 'FontSize', fontSize);

    end

    save(fullfile(savePath, 'superRes.mat'), 'TI', 'SI', 'EBSR', 'NI', 'BLI', 'interP');
    save(fullfile(savePath, 'EBSR.mat'), 'EBSR');

end
end


function cleanUpFun(param)

 saveData(param);
end