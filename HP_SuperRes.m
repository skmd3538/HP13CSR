function [ SR ] = HP_SuperRes( params )
%HP_SuperRes
%   Superresolution of hyperpolarized HP MR images
% The input params is a struct with the following fields:

% alpha controls tradeoff between matching low-resolution patch size
% and matching high-resolution neighbors (default:0.5444) 

% LRPatchSiz: low resolution patch size (default: 7)

% HRPatchSiz: high resolution patch size (default: 5)

% SRMethod: SuperRes algorithm (default: 1)

% externalLearning: use low resolution images from scanner (default: 0)

% compareToGroundTruth : compute PSNR and SSIM for sinc, NN and bicubic
% interpolations (default: 0)

% saveToDicom : save result to dicom format
% interpMethod: interpolation method. Same as defined for imresize
% Example usage:
% params.saveToDicom=1;params.SRMethod=1;params.externalLearning=1;
% params.compareToGroundTruth=1; HP_SuperRes(param)



if(~exist('params', 'var'))
    params = struct();
end

if(~isfield(params, 'externalLearning'))
    params.externalLearning = true;
end

if(~isfield(params, 'alpha'))
    params.alpha = 0.5444;
end

if(~isfield(params, 'interpMethod'))
    params.interpMethod = 'bicubic';
end

if(~isfield(params, 'rescale'))
    params.rescale = false;
end

if(~isfield(params, 'scale'))
    params.scale = 2;
%     params.alpha = 0.1923;
end

if(~isfield(params, 'bucketSize'))
    params.bucketSize = 50;
end

if(isfield(params, 'wtsFile'))
    wtsFile = params.wtsFile;
end

if(isfield(params, 'LRFile'))
    LRFile = params.LRFile;
    LRPath = params.LRPath;
end

if(~isfield(params, 'SRMethod'))
    params.SRMethod = 1;
end

if(~isfield(params, 'refImage'))
    params.refImage = 0;
end

str = {'Build Training Set', 'SuperRes'};
%datstr1 = {'training'};

if(~isfield(params, 'LRFile'))
    [selection,v] = listdlg('PromptString','Select data type:',...
        'SelectionMode','single',...
        'ListString',str,'ListSize',[160 60]);
    if v == 0; return; end
else
    selection = 2;
end


switch selection
    case {1}
        fprintf('\nSelect high res training images\n\n');
        [file,path] = uigetfile('*',...
            'Select high res training images', ...
            'MultiSelect', 'on');
        if isequal(file,0)
            disp('User selected Cancel')
        else
            wtName = split(strip(path, filesep), filesep);
            resWtFile = fullfile(path, strcat(wtName{end}, '_wts'));
            if ~iscell(file)
                file = {file};
            end
            if(params.externalLearning)
                fprintf('\nSelect external low res training images\n\n');
                [extfile,extpath] = uigetfile('*',...
                    'Select external low res training images', ...
                    'MultiSelect', 'on');
                if(isequal(extfile, 0)),return; end
                params.externalFile = fullfile(extpath, extfile);
            end
            %load an anatomic image if specified
             if(params.refImage)
                fprintf('\nSelect Anatomic Image\n\n');
                [extfile,extpath] = uigetfile('*',...
                    'Select Anatomic training images', ...
                    'MultiSelect', 'on');
                if(isequal(extfile, 0)),return; end
                params.refImage = fullfile(extpath, extfile);
            end
            buildTrainingSet(fullfile(path, file), resWtFile, params);
        end
    case 2
        if ~exist('wtsFile','var')
            fprintf('\nLoad weights file\n\n');
            [wtsNm,wtPath] = uigetfile('*.mat', 'Load weights file');
            if isequal(wtsNm, 0); return; end
            wtsFile = fullfile(wtPath, wtsNm);
        end
        if ~exist('LRFile', 'var')
            fprintf('\nSelect low res image\n\n');
            [LRFile, LRPath] = uigetfile('*', 'Select low res image', 'MultiSelect', 'on');
            if isequal(LRFile, 0); return; end
        end
        
        if(ischar(LRFile))
            inputImageFileNames = fullfile(LRPath,LRFile);
        end
        
        if(isfield(params, 'compareToGroundTruth') && params.compareToGroundTruth )
            fprintf('\nSelect ground truth image\n\n');
            [TruthFile, TruthPath] = uigetfile('*', 'Select ground truth image', 'MultiSelect', 'on');
            if ~isequal(TruthFile, 0), params.TruthImage = fullfile(TruthPath,TruthFile); end
        end
        
        if(isfield(params, 'saveToDicom') && params.saveToDicom)
            fprintf('\nSelect dicom template\n\n');
%             [DicomTemplate, DicomPath] = uigetfile('*', 'Select dicom template', 'MultiSelect', 'on');
            selpath = uigetdir('.','Select folder with dicom template');
            if ~isequal(selpath, 0), params.DicomTemplate = selpath; end
            %                 if(iscell(LRFile)), c = length(inputImageFileNames); else, c = 1; end
        end
        
        
        if(iscell(LRFile)), c = length(inputImageFileNames); else, c = 1; end
        for i = 1:c
            
            if(iscell(LRFile))
                inFile = inputImageFileNames{i}; 
            elseif(isnumeric(LRFile))
                inFile = LRFile;
            else
                inFile = inputImageFileNames; 
            end
            
            if(params.SRMethod == 1)
                if(~isnumeric(LRFile))
                    [~, fileNm, ~] = fileparts(LRFile);
                    savePath = fullfile(LRPath, sprintf('%s', fileNm));
                    params.savePath = savePath;
                    params.LRName = {fileNm};
                end
                [~, ~, superResImage] = ...
                    superRes(wtsFile, inFile, params);
            else
                params.inputImage = inFile;
                params.wtsFile = wtsFile;
                params.saveFile = LRPath;
                [superResImage] = MainDic_SR(params);
                params.MainDicImg = superResImage;
                %                         [superResImage] = ResidualDic_SR(params);
            end
            if(i == 1), SR = superResImage; else, SR = cat(4, SR, superResImage);  end
        end
        
end
        
end


