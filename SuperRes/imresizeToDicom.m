function imresizeToDicom( inputImageName, method, dicomTemplate )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
inputImage= load(inputImageName);
inputImage = abs(struct2array(inputImage(1)));

dic_info = dicominfo(dicomTemplate);
tgtSiz = [dic_info.Width dic_info.Height];


% dic_info.SeriesDescription = 'DHAbicubic';
dic_info.SeriesDescription = inputImageName;
dic_info.SeriesNumber = dic_info.SeriesNumber + 1;

resizedImg = imresize(inputImage, tgtSiz, method);

dicomwrite(mat2gray(rot90(resizedImg,1)), ...
    sprintf('%s_%s%d.dcm', inputImageName, method, tgtSiz(1)), dic_info)

end

