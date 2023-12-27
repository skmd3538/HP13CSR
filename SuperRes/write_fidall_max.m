function [ ] = write_fidall_max( )

[file,path] = uigetfile('*.mat');

msk_header.description = file;
msk_header.date = date;

if ~isequal(file,0)
    try
        load(fullfile(path, file), 'spectra_sum_fit', 'msk_header')
        spectra = spectra_sum_fit;
    catch
        load(fullfile(path, file), 'spec')
        spectra = spec;
    end
    max_image = shiftdim(max(abs(spectra),[],1),1);
    nm = [msk_header.description '_' strrep(msk_header.date, '/', '_')];
    nm(regexp(nm, '[\s,,/]'))=[];
    imwrite(mat2gray(max_image), fullfile(path,  strcat(nm, '_max_image.png')));
    dicomwrite(mat2gray(max_image), fullfile(path,  strcat(nm, '_max_image.dcm')));
end



end
