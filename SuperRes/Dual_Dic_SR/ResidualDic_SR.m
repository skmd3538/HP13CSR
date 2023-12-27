% function [ResSR] = ResidualDic_SR(InputImg, MainDicImg,ResidualDicImg, ResidualDic_MatName)
function [ResSR] = ResidualDic_SR(params)

B = [];
if(~exist('params', 'var'))
    params = struct();
end


if(isfield(params, 'MainDicImg'))
    MainIm = params.MainDicImg;
else
    %prompt for input file
end

if(isfield(params, 'saveFile'))
    savePath = params.saveFile;
else
    %prompt for weights file
end

% if(isfield(params, 'ResidualDic_MatName'))
%     ResidualDic_MatName = params.ResidualDic_MatName;
% else
%     %prompt for weights file
% end
[p, f, ~] = fileparts(params.wtsFile);
ResidualDic_MatName = sprintf('%s_ResDic.mat', fullfile(p,f));

load(ResidualDic_MatName);







%============================================

n=9; % block size
m=500; % number of atoms in the dictionary
dd=3; % margins in the image to avoid (dd*s to each side)
L=3; % number of atoms to use in the representation       


nd = ndims(MainIm);
sz = size(MainIm);
if nd > 2
    Yh = reshape(MainIm, sz(1), sz(2), []);
    sz3 = size(Yh,  3);
end

sz3 = size(Yh,  3);
ResSR = zeros(sz(1), sz(2), sz3);

        
for i = 1:sz3
    
    fprintf('Processing image frame:%d\n', i);
    Yhs = Yh(:,:,i);

Yl = im2double( Yhs );
 [N1, N2] = size(Yl);
% Extracting features
Yl1=conv2(Yl,[1,0,-1],'same'); % the filter is centered and scaled well for s=3
Yl2=conv2(Yl,[1,0,-1]','same');
Yl3=conv2(Yl,[1,0,-2,0,1]/2,'same');
Yl4=conv2(Yl,[1,0,-2,0,1]'/2,'same');

% Gathering the patches
clear Ptilde_l
Ptilde_l=zeros(4*n^2,(N1/s-2*dd)*(N2/s-2*dd)); 
counter=1;
for k1=s*dd+1:s:N1-s*dd
    for k2=s*dd+1:s:N2-s*dd
        Ptilde_l(:,counter)=[reshape(Yl1(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1]); ...
                                    reshape(Yl2(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1]);...
                                    reshape(Yl3(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1]);...
                                    reshape(Yl4(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1])];
        counter=counter+1;
    end;
end;

% Dimentionalily reduction
Pl=B'*Ptilde_l; 

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
    end;
end;
Yout=Yout./(Weight+1e-5)+Yl; 
Yout=min(max(Yout,0),255);

ResSR(:,:,i) = Yout;


end

[p, n, ~] = fileparts(params.inputImage);

savePath = fullfile(params.saveFile, sprintf('%s_DualDic_%dto%d', n, N1/s, N1));

ResSR = reshape(ResSR, sz);
savePath = sprintf('%s', fullfile(savePath, 'ResSR'));
save(savePath, 'ResSR');


% Yh = imresize(Yh, 2, 'nearest');
% Yout_name = strcat(ResIm,'_PSNR_',num2str(csnr(Yh,Yout,5,5)),'.png');
% imwrite(uint8(Yout),Yout_name);

end