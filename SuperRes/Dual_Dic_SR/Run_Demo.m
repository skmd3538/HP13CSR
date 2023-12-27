%   These MATLAB programs implement the Image Super-Resolution Algorithm via 
%   Dual Dictionary Learnning and Sparse Represenation as described in paper:
%   
%     J. Zhang, C. Zhao, R. Xiong, S. Ma, D. Zhao "Image Super-Resolution via 
%     Dual-Dictionary  Learning and Sparse Representation" in IEEE International 
%     Symposium of Circuits and Systems (ISCAS2012),pp. 1688¨C1691, Seoul, Korea, May 2012. 
%
% 
% -------------------------------------------------------------------------------------------------------
% The software implemented by MatLab 7.10.0(2010a) are included in this package.
%
% ------------------------------------------------------------------
% Requirements
% ------------------------------------------------------------------
% *) Matlab 7.10.0(2010a) or later with installed:
% ------------------------------------------------------------------
% Version 1.0
% Author: Jian Zhang
% Email:  jzhangcs@hit.edu.cn
% Last modified by J. Zhang, Feb. 2014

clear
clc

%Main_Dic_learning 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DicTrainingImage = 'HR-049.tif';
MainDic_MatName = 'Main_Dic.mat';
MainDic_ImName = 'HR-049_MainDic_SR.tif';
if ~exist('Main_Dic.mat') %#ok<EXIST>
    Main_Dic_Learning(DicTrainingImage,MainDic_MatName,MainDic_ImName);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Residual_Dic_learning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ResidualDic_MatName = 'Residual_Dic.mat';
if ~exist('Residual_Dic.mat') %#ok<EXIST> %
    Residual_Dic_Learning(DicTrainingImage,MainDic_ImName,ResidualDic_MatName);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ImgNo = 1:4
    
    switch ImgNo
        case 1
            InputImg = 'peppers.tif';
        case 2
            InputImg = 'Lena.tif';
        case 3
            InputImg = 'monarch512.tif';
        case 4
            InputImg = 'HR-067.tif';
    end
    
OrgX = double(imread(InputImg));

MainDicImg = strcat(InputImg,'_MainDic_SR');
MainSR = MainDic_SR(InputImg, MainDicImg,MainDic_MatName);
fprintf('PSNR of Super Resolution by Main Dictionary is %0.2f \n',csnr(OrgX,MainSR,5,5));

ResidualDicImg = strcat(InputImg,'_ResidualDic_SR');
ResSR = ResidualDic_SR(InputImg, MainDicImg,ResidualDicImg, ResidualDic_MatName);
fprintf('PSNR of Super Resolution by Main+Residual Dictionary is %0.2f \n',csnr(OrgX,ResSR,5,5));

end





%MSKCC test
DicTrainingImage = 'G32_M7.png';
MainDic_MatName = 'Main_Dic.mat';
MainDic_ImName = 'G32_M7_MainDic_SR.png';
if ~exist('Main_Dic.mat') %#ok<EXIST>
    Yout = Main_Dic_Learning(DicTrainingImage,MainDic_MatName,MainDic_ImName);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Residual_Dic_learning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ResidualDic_MatName = 'Residual_Dic.mat';
if ~exist('Residual_Dic.mat') %#ok<EXIST> %
    Residual_Dic_Learning(DicTrainingImage,MainDic_ImName,ResidualDic_MatName);
end

InputImg = 'G16_M7.png';

    
OrgX = double(imread(InputImg));

MainDicImg = strcat(InputImg,'_MainDic_SR');
MainSR = MainDic_SR(InputImg, MainDicImg,MainDic_MatName);
fprintf('PSNR of Super Resolution by Main Dictionary is %0.2f \n',csnr(OrgX,MainSR,5,5));

ResidualDicImg = strcat(InputImg,'_ResidualDic_SR');
ResSR = ResidualDic_SR(InputImg, MainDicImg,ResidualDicImg, ResidualDic_MatName);
fprintf('PSNR of Super Resolution by Main+Residual Dictionary is %0.2f \n',csnr(OrgX,ResSR,5,5));

imagesc(cat(2, G32, mat2gray(MainSR),  imresize(G16, 2, 'nearest'), imresize(G16, 2, 'bilinear'), imresize(G16, 2, 'bicubic')));daspect([1 1 1]);colorbar

MainSR = wiener2(MainSR,[10 10]);
MainSR = mat2gray(MainSR);
imagesc(cat(2, mat2gray(G32_M7),mat2gray(SI), mat2gray(NI), mat2gray(BCI), mat2gray(SPI), mat2gray(MainSR)));daspect([1 1 1]);colorbar;axis off
imagesc(cat(2, mat2gray(G64_M7),mat2gray(SI), mat2gray(NI), mat2gray(BCI), mat2gray(SPI), mat2gray(MainSR)));daspect([1 1 1]);colorbar;axis off
imagesc(cat(2, mat2gray(G128_M7),mat2gray(SI), mat2gray(NI), mat2gray(BCI), mat2gray(SPI), mat2gray(MainSR)));daspect([1 1 1]);colorbar;axis off

text(3, 31, 'Truth        SI         NI        BCI         SPI       SR');
T = table(mseSI, mseNI, mseBCI, mseSPI, mseSR);

imagesc(cat(2, mat2gray(G32_M12),mat2gray(MainSI), mat2gray(MainSR), mat2gray(MainBCI), mat2gray(MainSPI)));daspect([1 1 1]);colorbar;axis off

