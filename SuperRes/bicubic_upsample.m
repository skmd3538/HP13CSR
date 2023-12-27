inputImageName = 'MetMaps_VitaminC_MAG_INTEGRATEI1.dcm';
inputImage = dicomread(inputImageName);
dic_info = dicominfo(inputImageName);
tgtSiz = [48 48];

dic_info.PixelSpacing(1) = dic_info.PixelSpacing(1) * double(dic_info.Width)/tgtSiz(1);
dic_info.PixelSpacing(2) = dic_info.PixelSpacing(2) * double(dic_info.Height)/tgtSiz(2);
dic_info.Width = tgtSiz(1);
dic_info.Height = tgtSiz(2);
dic_info.Rows = tgtSiz(1);
dic_info.Columns = tgtSiz(2);
% dic_info.SeriesDescription = 'DHAbicubic';
dic_info.SeriesDescription = 'VitCbicubic';
dic_info.SeriesNumber = dic_info.SeriesNumber + 1;

Bicubic48 = imresize(inputImage, [48 48], 'bicubic');
dicomwrite(mat2gray(rot90(Bicubic48,1)), sprintf('M1195_zf_VitCbicubic.dcm'), dic_info)
% dicomwrite(mat2gray(rot90(Bicubic48,1)), sprintf('M1195_DHAbicubic.dcm'), dic_info)
% nm = split(inputImageName, '.');
% nm = nm{1};
% dicomwrite(mat2gray(rot90(Bicubic48,1)), sprintf('M1195_DHAbicubic.dcm'), dic_info)
% dicomwrite(mat2gray(Bicubic48), sprintf('M1195_VitCbicubic.dcm'), dic_info)



EBSR = dicomread('M1195_VitaminC_SR24_SR48.dcm');
    dicomwrite(mat2gray(rot90(EBSR,1)), 'M1195_VitaminC_SR24_SR48.dcm');
