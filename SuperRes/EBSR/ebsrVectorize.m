%performs concatenation of 49 low res pixels to 9 high res pixels
%for more info, see top comments in createVectors.m
function output = ebsrVectorize(alpha, lowResPatches, highResPatches, i, j)

lowResPatch = accessPatch(lowResPatches, i,j);  %7x7x3
highResPatch = accessPatch(highResPatches, i+1,j+1);

hRPatchLen = size(highResPatches, 3);
lRPatchLen = size(lowResPatches, 3);
patchLen = (lRPatchLen * lRPatchLen + (hRPatchLen + hRPatchLen - 1));

output = zeros(1,patchLen,3);

for k = 1:3
    highResOutput = [];
    highResOutput = [highResPatch(1,:,k) highResPatch(2:end,1,k)'];
    highResOutput = highResOutput * alpha;
    output(:,:,k) = [reshape(lowResPatch(:,:,k), 1, []) highResOutput];
end

output = reshape(output, 1, 3*patchLen);

end

