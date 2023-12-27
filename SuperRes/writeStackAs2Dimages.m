metabolite_maps_fit = permute(metabolite_maps_fit, [2 3 1]);
StackImage = metabolite_maps_fit;
% mkdir TrainingImages
% cd TrainingImages\
StackImage = CorT2S_EPI20x40x20_5;
for i=1:size(StackImage,3)
    X = mat2gray(StackImage(:,:,i));
     save(sprintf('TrainImage%d',i),'X');
%     imwrite(double(X), sprintf('TrainImage%d.tif',i));
end