function [Yout]=Main_Dic_Learning(params)

%This function is used to learn the main dictionary 
%(including high and low dictionary) associated with
%different image, and save them 
%2010.11.5
%============================================
%                              Training the two dictionaries
%============================================
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

if(isfield(params, 'wtsFile'))
    MatName = params.wtsFile;
else
    %prompt for weights file
end


% Setting parameters
n=9; % block size
%m=500; % number of atoms in the dictionary
m=100; % number of atoms in the dictionary

dd=3; % margins in the image to avoid (dd*s to each side)
L=3; % number of atoms to use in the representation

savefile = MatName;

% Preparing the low-and high resolution images
[~, ~, fExt] = fileparts(ImageIN);
if(isdicom(ImageIN))
    Yh = im2double(dicomread(ImageIN));
elseif strcmp(fExt, '.mat')
    Yh = load(ImageIN);
    Yh = abs(struct2array(Yh(1)));
else
    Yh = imread(ImageIN);
end


[N1, N2] = size(Yh);
m = N2;
N1=floor(N1/s)*s; 
N2=floor(N2/s)*s;
Yh=im2double(Yh(1:N1,1:N2)); % so that it dcales down to an integer size
Yh=Yh*255; 

% Creating the low-resolution image to train on 
% Using different blur operators
%Zl=conv2(Yh,[1 2 1]/4,'same');
%Zl=conv2(Zl,[1 2 1]'/4,'same');

HH = fspecial('Gaussian',[5 5],1);
Zl = imfilter(Yh,HH,'same');
Zl=Zl(1:s:end,1:s:end); 

%Zl = imresize(Yh,1/2,'bicubic');


% Upscaling Zl to the original resolution
%N=size(Yh,1);
[posY,posX]=meshgrid(1:s:N2,1:s:N1); 
[posY0,posX0]=meshgrid(1:1:N2,1:1:N1); 
Yl=interp2(posY,posX,Zl,posY0,posX0,'bicubic');

%Yl = imresize(Zl,2,'bicubic');

% Extracting features
Eh=Yh-Yl; % pre-processing of the high-resolution image
Yl1=conv2(Yl,[1,0,-1],'same'); % the filter is centered and scaled well for s=3
Yl2=conv2(Yl,[1,0,-1]','same');
Yl3=conv2(Yl,[1,0,-2,0,1]/2,'same');
Yl4=conv2(Yl,[1,0,-2,0,1]'/2,'same');

% Gathering the patches
Ph=zeros(n^2,(N1/s-2*dd)*(N2/s-2*dd)); 
Ptilde_l=zeros(4*n^2,(N1/s-2*dd)*(N2/s-2*dd)); 
n2=(n-1)/2; 
counter=1;
for k1=s*dd+1:s:N1-s*dd
    for k2=s*dd+1:s:N2-s*dd
        Ph(:,counter)=reshape(Eh(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1]);
        
         Ptilde_l(:,counter)=[reshape(Yl1(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1]); ...
                                    reshape(Yl2(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1]);...
                                    reshape(Yl3(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1]);...
                                    reshape(Yl4(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1])];
        counter=counter+1;
    end
end

% Dimentionalily reduction
R=Ptilde_l*Ptilde_l'; 
[B,SS]=eig(R); 
Permute=fliplr(eye(size(R,1))); 
SS=Permute*SS*Permute; % so that the eigenvalues are sorted descending
B=B*Permute; 
energy=cumsum(diag(SS))/sum(diag(SS)); 
% figure(1); clf; plot(energy)
pos=find(energy>0.999,1);
B=B(:,1:pos);
disp(['Effective dimension: ',num2str(pos)]); 
disp(['The relative error is: ',...
        num2str(mean(sum((B*B'*Ptilde_l-Ptilde_l).^2))/...
        mean(sum((Ptilde_l).^2)))]); % showing the relative error
Pl=B'*Ptilde_l; 

% Low-Res. Dictionary Learning 
param.errorFlag=0;
param.K=m; 
%param.K=size(Pl, 2);
param.numIteration=20; 
param.InitializationMethod='DataElements'; 
param.TrueDictionary=0;
param.Method='KSVD';
param.L=L;   
[Al,output]=TrainDic_Fast(Pl,param);
Q=output.CoefMatrix; 

% High-Resolution Dictionary Leanring
Ah=Ph*Q'*inv(Q*Q');

%============================================
%                  Sanity Check - Interpolating the training image
%============================================

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
%csnr(Yh,Yl,5,5)
%csnr(Yh,Yout,5,5)

savefile = sprintf('%s_%dto%d.mat', savefile, N1/s, N1);
save(savefile,'Al','Ah','B','dd','L','n','n2','m','s');
fprintf('Weights saved to %s\n', savefile);
% imwrite(uint8(Yout),Main_dic_ImName,'tif');

end