function buildTrainingSet(inputImageFileNames, saveFileName, params)

%c contains number of input files
[~, c] = size(inputImageFileNames);

%training set is a huge set of key-value pairs
keys = [];
values = [];


if params.SRMethod == 1
for i = 1:c
    display(['Generating training set from image:  ' inputImageFileNames{i}]);
    if(isfield(params, 'externalFile'))
        if c == 1, params.externalFileName = params.externalFile; else, params.externalFileName = params.externalFile{i}; end
        display(['Generating training set with external file:  ' params.externalFileName]);
    end

    
    [tKeys, tValues, imgSiz, tgtSiz, params] = train(inputImageFileNames{i}, params);
   
    keys = [keys ; tKeys];
    values = [values ; tValues];
end
LRPatchSiz = params.LRPatchSiz;
HRPatchSiz = params.HRPatchSiz;
[~, fileNm, ~] = fileparts(inputImageFileNames{1});
alphaStr = strrep(sprintf('%.3f', params.alpha), '.', 'P');
alphaStr = sprintf('%sLRP%dHRP%d', alphaStr, params.LRPatchSiz, params.HRPatchSiz);
saveFileName = sprintf('%s_%s_EBSR_alpha%s_%dx%dto%dx%d.mat', saveFileName,  fileNm, alphaStr, imgSiz(1),imgSiz(2), tgtSiz(1), tgtSiz(2));
save(saveFileName, 'keys', 'values', 'LRPatchSiz', 'HRPatchSiz');
fprintf('Weights saved to %s\n', saveFileName);
elseif params.SRMethod == 2
    inputImageFileNames = inputImageFileNames{1};
    [~, n, ~] = fileparts( inputImageFileNames);
    saveFileName = sprintf('%s_DualDic_%s', saveFileName, n);
    params.wtsFile = saveFileName;
    params.inputImage = inputImageFileNames;
    Yout = Main_Dic_Learning(params);
     params.Main_Dic_Im = Yout;
    Residual_Dic_Learning(params);
end

