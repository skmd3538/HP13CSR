
function [MainSR] = MainDic_SR(params)

printf(' Super Resolution Using Main Dictionary:%s\n', params.wtsFile);


if(~exist('params', 'var'))
    params = struct();
end

% scale-down factor
if(isfield(params, 'scale'))
    s = params.scale;
else
   s = 2; 
end

if(isfield(params, 'inputImage'))
    ImageIN = params.inputImage;
else
    %prompt for input file
end

if(isfield(params, 'saveFile'))
    saveFile = params.saveFile;
else
    %prompt for weights file
end


if(isfield(params, 'wtsFile'))
    MainDic_MatName = params.wtsFile;
else
    %prompt for weights file
end

if(isfield(params, 'TruthImage'))
    TruthImage = params.TruthImage;
else
    %prompt for weights file
end



B = [];
load(MainDic_MatName);

fprintf('Input Image: %s\n',ImageIN);
%============================================

% Setting parameters
n=9; % block size

m= 500; % number of atoms in the dictionary
dd=3; % margins in the image to avoid (dd*s to each side)
L=3; % number of atoms to use in the representation
n2=(n-1)/2;

[~, ~, fExt] = fileparts(ImageIN);
if(isdicom(ImageIN))
    Yh = im2double(dicomread(ImageIN));
elseif strcmp(fExt, '.mat')
    Yh = load(ImageIN);
    Yh = abs(struct2array(Yh(1)));
nd = ndims(Yh);
        sz = size(Yh);
    if nd > 2
        if(isfield(params, 'tp'))
             Yh = reshape(Yh, [], sz(end));
             Yh = Yh(:, params.tp);
             sz(end) = numel(params.tp);
        end
        Yh = reshape(Yh, sz(1), sz(2), []);
    end
else
    Yh = imread(ImageIN);
end

sz(1) = size(Yh, 1);
sz(2) = size(Yh, 2);
sz3 = size(Yh,  3);

% Yh = imresize(Yh, 2, 'nearest');
MainSR = zeros(sz(1)*s, sz(2)*s, sz3);
MainSI = zeros(sz(1)*s, sz(2)*s, sz3);
MainBCI = zeros(sz(1)*s, sz(2)*s, sz3);
MainBLI = zeros(sz(1)*s, sz(2)*s, sz3);
MainSPI = zeros(sz(1)*s, sz(2)*s, sz3);
MainNI = zeros(sz(1)*s, sz(2)*s, sz3);
        
for i = 1:sz3
    
    fprintf('Processing image frame:%d\n', i);
    Yhs = Yh(:,:,i);
    [N1, N2] = size(Yhs);
    N1=floor(N1/s)*s;
    N2=floor(N2/s)*s;
    Yhs=im2double(Yhs(1:N1,1:N2)); % so that it dcales down to an integer size
    Yhs=Yhs*255;
    Zl = Yhs;
    N1 = N1 * s;
    N2 = N2 * s;
    
    % Upscaling Zl to the original resolution
    [posY,posX]=meshgrid(1:s:N2,1:s:N1);
    [posY0,posX0]=meshgrid(1:1:N2,1:1:N1);
    Yl=interp2(posY,posX,Zl,posY0,posX0,'bicubic');
    SPI=interp2(posY,posX,Zl,posY0,posX0,'spline');
    NI=interp2(posY,posX,Zl,posY0,posX0,'nearest');
    BLI=interp2(posY,posX,Zl,posY0,posX0,'bilinear');
    SI = sincinterp (Yhs, N1, N1); %sinc interpolation
    BCI = Yl;

    
    
    % Extracting features
    Yl1=conv2(Yl,[1,0,-1],'same'); % the filter is centered and scaled well for s=3
    Yl2=conv2(Yl,[1,0,-1]','same');
    Yl3=conv2(Yl,[1,0,-2,0,1]/2,'same');
    Yl4=conv2(Yl,[1,0,-2,0,1]'/2,'same');
    
    % Gathering the patches
    Ptilde_l=zeros(4*n^2,(N1/s-2*dd)*(N2/s-2*dd));
    counter=1;
    for k1=s*dd+1:s:N1-s*dd
        for k2=s*dd+1:s:N2-s*dd
            Ptilde_l(:,counter)=[reshape(Yl1(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1]); ...
                reshape(Yl2(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1]);...
                reshape(Yl3(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1]);...
                reshape(Yl4(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1])];
            counter=counter+1;
        end
    end
    
    % Dimentionalily reduction
    Bvar = B;
    Pl=Bvar'*Ptilde_l;
    
    % Cleaning up some space in memory
    clear posX posY posX0 posY0 Zl Yi4 Yl3 Yl2 Yl1 Ptilde_l
    
    % Sparse coding of the low-res patches
    % Does not run because of memory problems
    Pl_size = size(Pl,2);
    if Pl_size<20000
        Q=omp(Al'*Pl,Al'*Al,L);
    elseif Pl_size<40000
        Q1=omp(Al'*Pl(:,1:20000),Al'*Al,L);
        Q2=omp(Al'*Pl(:,20001:end),Al'*Al,L);
        Q=[Q1,Q2];
    elseif Pl_size<60000
        Q1=omp(Al'*Pl(:,1:20000),Al'*Al,L);
        Q2=omp(Al'*Pl(:,20001:40000),Al'*Al,L);
        Q3=omp(Al'*Pl(:,40001:end),Al'*Al,L);
        Q=[Q1,Q2,Q3];
    elseif Pl_size<80000
        Q1=omp(Al'*Pl(:,1:20000),Al'*Al,L);
        Q2=omp(Al'*Pl(:,20001:40000),Al'*Al,L);
        Q3=omp(Al'*Pl(:,40001:60000),Al'*Al,L);
        Q4=omp(Al'*Pl(:,60001:end),Al'*Al,L);
        Q=[Q1,Q2,Q3,Q4];
    elseif Pl_size<100000
        Q1=omp(Al'*Pl(:,1:20000),Al'*Al,L);
        Q2=omp(Al'*Pl(:,20001:40000),Al'*Al,L);
        Q3=omp(Al'*Pl(:,40001:60000),Al'*Al,L);
        Q4=omp(Al'*Pl(:,60001:80000),Al'*Al,L);
        Q5=omp(Al'*Pl(:,80001:end),Al'*Al,L);
        Q=[Q1,Q2,Q3,Q4,Q5];
    end
    
    
    % Recover the image
    Ph_hat=Ah*Q;
    Yout=Yl*0;
    Weight=Yl*0;
    counter=1;
    for k1=s*dd+1:s:N1-s*dd
        for k2=s*dd+1:s:N2-s*dd
            patch=reshape(Ph_hat(:,counter),[n,n]);
             Yout(k1-n2:k1+n2,k2-n2:k2+n2)=...
                Yout(k1-n2:k1+n2,k2-n2:k2+n2)+patch;
            Weight(k1-n2:k1+n2,k2-n2:k2+n2)=...
                Weight(k1-n2:k1+n2,k2-n2:k2+n2)+1;
            counter=counter+1;
        end
    end
    Yout=Yout./(Weight+1e-5)+Yl;
    Yout=min(max(Yout,0),255);
    
    
    MainSR(:,:,i) = Yout;
    MainSI(:,:,i) = SI;
    MainBCI(:,:,i) = BCI;
    MainBLI(:,:,i) = BLI;
    MainSPI(:,:,i) = SPI;
    MainNI(:,:,i) = NI;

    
end
% Yh = imresize(Yh, 2, 'nearest');
% 
% Yl_name = strcat(ImageIN,'_bicubic_PSNR_',num2str(csnr(Yh,Yl,5,5)),'.png');
% Yout_name = strcat(MainIm,'_PSNR_',num2str(csnr(Yh,Yout,5,5)),'.png');
% imwrite(uint8(Yl),Yl_name);
% imwrite(uint8(Yout),MainIm,'tif');
% imwrite(uint8(Yout),Yout_name);

[~, b, ~] = fileparts(ImageIN);
savePath = fullfile(saveFile, sprintf('%s_DualDic_%dto%d', b, N1/s, N1));
if ~exist(savePath, 'dir')
    mkdir(savePath)
end

sz(1) = sz(1) * s;
sz(2) = sz(2) * s;
MainSR = reshape(MainSR, sz);
MainSI = reshape(MainSI, sz);
MainBCI = reshape(MainBCI, sz);
MainBLI = reshape(MainBLI, sz);
MainSPI = reshape(MainSPI, sz);
MainNI = reshape(MainNI, sz);

save(fullfile(savePath, 'superRes.mat'), ...
    'MainSR', 'MainSI','MainBCI', 'MainBLI', 'MainSPI', 'MainNI');

params.saveFile = savePath;

if exist('TruthImage', 'var')  
    
    
    [~, ~, fExt] = fileparts(TruthImage);
    if(isdicom(TruthImage))
        TI  = im2double(dicomread(TruthImage));
    elseif strcmp(fExt, '.mat')
        TI  = load(TruthImage);
        TI  = im2double(abs(struct2array(TI (1))));
    else
        TI = im2double(imread(TruthImage));
    end
    
    NI(isnan(NI)) = 0;
    BCI(isnan(BCI)) = 0;
    
    mseSR = immse(TI,MainSR);
    mseSI = immse(TI,SI);
    mseBCI = immse(TI,BCI);
    mseBLI = immse(TI,BLI);
    mseSPI = immse(TI,SPI);
    mseNI = immse(TI,NI);
    
    ssimSR = ssim(TI,MainSR);
    ssimSI = ssim(TI,SI);
    ssimBCI = ssim(TI,BCI);
    ssimSPI = ssim(TI,SPI);
    ssimNI = ssim(TI,NI);
    
    cmap = colormap(jet(64));
    fontSize = 20;
    
    subplot(1,5,1),imshow(mat2gray(TI), 'DisplayRange', [],'ColorMap', cmap),title('Truth','FontSize', fontSize);hold  on
    subplot(1,5,2),imshow(mat2gray(SI), 'DisplayRange', [],'ColorMap', cmap),title('Zerofill','FontSize', fontSize);
    xlabel(sprintf('MSE:%.01f, SSIM:%.03f', mseSI, ssimSI), 'FontSize', fontSize);
    subplot(1,5,3),imshow(mat2gray(NI), 'DisplayRange', [],'ColorMap', cmap),title('Nearest','FontSize', fontSize);
    xlabel(sprintf('MSE:%.01f, SSIM:%.03f', mseNI, ssimNI), 'FontSize', fontSize);
    subplot(1,5,4),imshow(mat2gray(BCI), 'DisplayRange', [],'ColorMap', cmap),title('Bicubic','FontSize', fontSize);
    xlabel(sprintf('MSE:%.01f, SSIM:%.03f', mseBCI, ssimBCI), 'FontSize', fontSize);
    subplot(1,5,5),imshow(mat2gray(MainSR), 'DisplayRange', [],'ColorMap', cmap),title('DualDicRes','FontSize', fontSize);hold off
    xlabel(sprintf('MSE:%.01f, SSIM:%.03f', mseSR, ssimSR),  'FontSize', fontSize);
    
    
    
    
    
%     
%     
%     imagesc(cat(2, mat2gray(TI),mat2gray(NI), mat2gray(SI), ...
%         mat2gray(BCI), mat2gray(BLI), mat2gray(SPI), mat2gray(MainSR)));
%     daspect([1 1 1]);colorbar;axis off
%     title(sprintf('Truth    Nearest     Sinc    Bicubic     Bilinear    Spline      SuperRes'));
    
    save(fullfile(savePath, 'superRes.mat'), ...
    'mseSR', 'mseSI','mseBCI', 'mseBLI', 'mseSPI', 'mseNI', ...
    'ssimSR', 'ssimSI','ssimBCI', 'ssimSPI', 'ssimNI', '-append');
    
end


end