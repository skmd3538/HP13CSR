% function []=Residual_Dic_Learning(ImageIN,Main_dic_ImName,MatName)
function []=Residual_Dic_Learning(params)

%This function is used to learn the residual dictionary 
%(including high and low dictionary) associated with
%different image, and save them 
%2010.11.5

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

if(isfield(params, 'Main_Dic_Im'))
    Main_dic_Im = params.Main_Dic_Im;
else
    %prompt for main dic. image
end

% Setting parameters
n=9; % block size
% m=500; % number of atoms in the dictionary
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

[N1 N2] = size(Yh);
N1=floor(N1/s)*s; 
N2=floor(N2/s)*s;
Yh=im2double(Yh(1:N1,1:N2)); % so that it dcales down to an integer size
Yh=Yh*255; 

%Yl = double(imread(Main_dic_ImName));
Yl = Main_dic_Im;
Yl=im2double(Yl(1:N1,1:N2)); 

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
    end;
end;

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
param.numIteration=20; 
param.InitializationMethod='DataElements'; 
param.TrueDictionary=0;
param.Method='KSVD';
param.L=L;   
[Al,output]=TrainDic_Fast(Pl,param);
Q=output.CoefMatrix; 

% High-Resolution Dictionary Leanring
Ah=Ph*Q'*inv(Q*Q');
savefile = sprintf('%s_%dto%d_ResDic.mat', savefile, N1/s, N1);
save(savefile,'Al','Ah','B','dd','L','n','n2','m','s');

end