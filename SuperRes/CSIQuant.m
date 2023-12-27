clear
close all
load('lb0_zf0-0-0-0.mat', 'msk_header','spectra_fit', 'fit_cs')
% load('lb20_zf0-0-1-0.mat')
% load('lb0_zf0-0-0-0.mat', 'msk_header', 'spectra_fit', 'spectra_bcpc', 'fit_cs')
% Substrate
% figure_name = 'M1195';
% name_substrate = 'DHA';
% name_met1 = 'VitaminC';
% name_met2 = 'Urea';
figure_name = 'HP0032';
name_substrate = 'Pyruvate';
name_met1 = 'Lactate';
name_met2 = 'Alanine';
% 
figure_name = 'Phantom';
name_substrate = 'Urea';
name_met1 = 'Alanine';
name_met2 = 'Acetate';

dim_spect = 1;
dim_read = 2;
dim_phase = 3;
dim_slice = 4;
dim_time = 5;
dim_echo = 6;
dim_coil = 7;
dim_ave = 8;
spectra_proc = spectra_bcpc;
label_disp = 'MIP spectrum (mag)'; spect_disp = max( max( max( max( max( max( max( abs( spectra_proc ),[],dim_read ),[],dim_phase ),[],dim_slice ),[],dim_time ),[],6 ),[],7 ),[],8 ); % MIP spectrum
[ disp_max,disp_idx ] = max( spect_disp );
% msk_header.center_freq = h.rdb_hdr.ps_mps_freq/1e7;

chemshift = msk_header2chemshift_axis( msk_header );




ppm_substrate = [fit_cs(1)-2 fit_cs(1)+2]
fit_cs_orig = fit_cs;
if(length(fit_cs) == 1)
    fit_cs(2) = 200;
    fit_cs(3) = 210;
end

ppm_met1 = [fit_cs(2)-2 fit_cs(2)+2]
ppm_met2 = [fit_cs(3)-1 fit_cs(3)+1]
% ppm_met2 = [162.1000 182.8000]

cmap = hot;
alpha_transparency = 0.4;
% ppm_substrate = [fit_cs(1)-1 fit_cs(3)+1]

data = spectra_proc;
idx_substrate = intersect( find( chemshift>=min( ppm_substrate )),find( chemshift<=max( ppm_substrate )));
idx_met1 = intersect( find( chemshift>=min( ppm_met1 )),find( chemshift<=max( ppm_met1 )));
idx_met2 = intersect( find( chemshift>=min( ppm_met2 )),find( chemshift<=max( ppm_met2 )));
ratio_map = squeeze( sum( abs( data( idx_met1,:,: )),1 )./sum( abs( data( idx_substrate,:,: )),1 ));
overlay_transparency = alpha_transparency*ones( size( ratio_map ));

spect_substrate = zeros(length(spect_disp));
spect_substrate(idx_substrate) = spect_disp(idx_substrate);

spect_met = zeros(length(spect_disp));
spect_met(idx_met1) = spect_disp(idx_met1);

spect_met2 = zeros(length(spect_disp));
spect_met2(idx_met2) = spect_disp(idx_met2);


h = plot_spectrum( chemshift,spect_disp ); title( label_disp ); hold on; xlabel( 'chemical shift [ppm]' );
h = plot_spectrum( chemshift,spect_substrate ); title( label_disp ); hold on; xlabel( 'chemical shift [ppm]' );
h = plot_spectrum( chemshift,spect_met ); title( label_disp ); hold on; xlabel( 'chemical shift [ppm]' );
h = plot_spectrum( chemshift,spect_met2 ); title( label_disp ); hold on; xlabel( 'chemical shift [ppm]' );
hold on


substrate_scale = 1;
map_substrate = squeeze( sum( abs( data( idx_substrate,:,: )),1 ));
map_met1 = squeeze( sum( abs( data( idx_met1,:,: )),1 ));
map_substrate_rgb = ind2rgb( round( map_substrate*length( cmap )/max( map_substrate(:))),cmap );
map_met1_rgb = ind2rgb( round(  map_met1*length( cmap )/max( map_substrate(:)/substrate_scale )),cmap );

map_met2 = squeeze( sum( abs( data( idx_met2,:,: )),1 ));
map_met2_rgb = ind2rgb( round(  map_met2*length( cmap )/max( map_substrate(:)/substrate_scale )),cmap );


pathname = pwd;
display_params.fig_label = [ figure_name,' - ',name_substrate ];
display_params.color_overlay = cat( 3,map_substrate_rgb,overlay_transparency );
disp_spectra_overlay( max( abs( data ),[],dim_time ),msk_header,pathname,display_params );
set( gcf,'Renderer','painters' ); saveas( gcf,fullfile( pathname,sprintf( '%s_%s-map_%0.0fppm',figure_name,name_substrate,mean( ppm_substrate ))),'svg' );

% Met1
display_params.fig_label = sprintf( '%s - %s x%0.2f',figure_name,name_met1,substrate_scale );
display_params.color_overlay = cat( 3,map_met1_rgb,overlay_transparency );
disp_spectra_overlay( max( abs( data ),[],dim_time ),msk_header,pathname,display_params );


% Met2
display_params.fig_label = sprintf( '%s - %s x%0.2f',figure_name,name_met2,substrate_scale );
display_params.color_overlay = cat( 3,map_met2_rgb,overlay_transparency );
disp_spectra_overlay( max( abs( data ),[],dim_time ),msk_header,pathname,display_params );

scan = split(pwd, '\');
scan = scan{end-1};
nm = sprintf('%s_%s', msk_header.description, scan);
nm = strrep(strrep(strrep(nm, '(', ''), ')', ''), ' ', '_');
map_substrate = mat2gray(map_substrate);
map_met1 = mat2gray(map_met1);
map_met2 = mat2gray(map_met2);

save(sprintf('%s_%s', nm, name_substrate), 'map_substrate');
copyfile(sprintf('%s_%s.mat', nm, name_substrate), '../../')
if(length(fit_cs_orig) > 1)
    save(sprintf('%s_%s', nm, name_met1), 'map_met1');
    copyfile(sprintf('%s_%s.mat', nm, name_met1), '../../')
    save(sprintf('%s_%s', nm, name_met2), 'map_met2');
    copyfile(sprintf('%s_%s.mat', nm, name_met2), '../../')
end

% close all
% clear