metMaps = abs(IDEALMaps.metaboliteMaps);
sz = size(metMaps);
metMaps = reshape(metMaps, [sz(1), sz(2), prod(sz(3:5))]);

for i = 1:size(metMaps, 3)
    metMaps(:,:,i) = mat2gray(imnlmfilt(metMaps(:,:,i)));
end

metMapsGray = reshape(metMaps, sz);
