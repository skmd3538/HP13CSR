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

substrate_ppm = 143;
phantom_ppm = 182;
figure; hold on;
h = plot_spectrum( chemshift,spect_disp ); title( label_disp ); hold on; xlabel( 'chemical shift [ppm]' );
set( h,'Color','k' );
 plot( chemshift( disp_idx ),spect_disp( disp_idx ),'ro' );
                text( chemshift( disp_idx )-1,double(disp_max),sprintf( '(%0.2f,%0.2g)',chemshift( disp_idx ),double(disp_max) ),'Color','r' );
cs_disp = sort([ substrate_ppm,phantom_ppm ],'ascend' )+[ -20,20 ];
                xlim([ max([ min( chemshift ),min( cs_disp )]),min([ max( chemshift ),max( cs_disp )])]);
peak_center = 146.2;


[ ~,idx_center ] = min( abs( peak_center-chemshift ));
ax1 = subplot( 2,1,1 ); h = plot_spectrum( chemshift,spect_disp ); title( label_disp ); xlabel( 'chemical shift [ppm]' ); hold on;
set( h,'Color','k' );
subplot( 2,1,1 ); plot( chemshift( idx_center ),spect_disp( idx_center ),'or' ); legend([{ 'unshifted spectrum' },{ 'substrate center' }],'Location','best','Box','off' );
 peak_ppm = 163;

if( isfield( msk_header,'matrix_size_recon' ))
    spect_shift = round(( peak_center-peak_ppm )/( msk_header.sweep_width/msk_header.center_freq/msk_header.matrix_size_recon( dim_spect )));
else
    spect_shift = round(( peak_center-peak_ppm )/( msk_header.sweep_width/msk_header.center_freq/msk_header.pts_spect_recon ));
end
spect_disp = circshift( spect_disp,spect_shift,1 );
%         spectra_8d = circshift( spectra_8d,spect_shift,1 );
% noisy_spectra_shift = circshift( noisy_spectra,spect_shift,1 );

ax2 = subplot( 2,1,2 ); h = plot_spectrum( chemshift,spect_disp ); legend({ 'shifted spectrum' },'Location','best','Box','off' );
set( h,'Color','k' );
% linkaxes([ ax1,ax2 ],'xy' ); xlim([ min( chemshift ),max( chemshift )]);
xlabel( 'chemical shift [ppm]' );
% set( gcf,'Name',fig_label )

% ppm_substrate = [169.5 172.5];
% ppm_met1 = [182.7 184.3];
% ppm_substrate = [150 190]; %phantom data
% ppm_met1 = [140 150];
% ppm_substrate = [172 176];
% ppm_met1 = [178.7 177.2];
% ppm_substrate = [169.5 172.5];
% ppm_met1 = [182.7 184.3];
% fit_cs = [-16.7700   -0.1189    5.5890]
% ppm_substrate = [fit_cs(1)-2 fit_cs(1)+2];
% ppm_met1 = [fit_cs(2)-2 fit_cs(2)+1];
% ppm_met2 = [fit_cs(3)-1 fit_cs(3)+1];

% ppm_substrate = [163-2  163+2];
% ppm_met1 = [178-2 178+3]
% ppm_met2 = [185-1 185+1]
ppm_substrate = [fit_cs(1)-2 fit_cs(1)+1.5];
ppm_met1 = [fit_cs(2)-1 fit_cs(2)+1.5]
ppm_met2 = [fit_cs(3)-2 fit_cs(3)+1.5]
    cmap = hot;
    alpha_transparency = 0.4;
    
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

% Substrate
figure_name = 'Phantom';
name_substrate = 'DHA';
name_met1 = 'Vitamin C';
name_met2 = 'Urea';

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

save M1075_20201116_DHA_CSI_PBS_1wk_DHA map_substrate
save M1075_20201116_DHA_CSI_PBS_1wk_VitC map_met1
save M1075_20201116_DHA_CSI_PBS_1wk_Urea map_met2

    
    
        subplot(1,6,1),imshow(mat2gray(fliplr(SI))),title('Zerofill'); colormap('jet')
    subplot(1,6,2),imshow(mat2gray(NI')),title('Nearest');colormap('jet')
    subplot(1,6,3),imshow(mat2gray(interP')),title('Bicubic');colormap('jet')
    subplot(1,6,4),imshow(mat2gray(EBSR')),title('EBSR');colormap('jet')
    subplot(1,6,5),imshow(mat2gray(Enmd')),title('EBSR(medfilt)');colormap('jet')
    subplot(1,6,6),imshow(mat2gray(Enlm')),title('EBSR(nlmfilt)');colormap('jet')