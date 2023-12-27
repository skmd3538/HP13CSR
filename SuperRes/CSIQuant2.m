load('lb0_zf0-0-0-0.mat', 'msk_header', 'spectra_bcpc', 'fit_cs')
spectra_proc = spectra_bcpc;
dim_spect = 1;dim_read = 2;dim_phase = 3;dim_slice = 4;dim_time = 5;dim_echo = 6;dim_coil = 7;dim_ave = 8;
label_disp = 'MIP spectrum (mag)'; spect_disp = max( max( max( max( max( max( max( abs( spectra_proc ),[],dim_read ),[],dim_phase ),[],dim_slice ),[],dim_time ),[],6 ),[],7 ),[],8 ); % MIP spectrum
[ disp_max,disp_idx ] = max( spect_disp );
chemshift = msk_header2chemshift_axis( msk_header );

        
figure; hold on;
h = plot_spectrum( chemshift,spect_disp ); title( label_disp ); hold on; xlabel( 'chemical shift [ppm]' );
set( h,'Color','k' );

% ppm_substrate = [169.5 172.5];
% ppm_met1 = [182.7 184.3];
cmap = hot;
alpha_transparency = 0.4;
    
ppm_substrate = [fit_cs(2)-2 fit_cs(2)+2]; %phantom data


   
    
data = spectra_bcpc;
idx_substrate = intersect( find( chemshift>=min( ppm_substrate )),find( chemshift<=max( ppm_substrate )));


spect_substrate = zeros(length(spect_disp));
spect_substrate(idx_substrate) = spect_disp(idx_substrate);


h = plot_spectrum( chemshift,spect_disp ); title( label_disp ); hold on; xlabel( 'chemical shift [ppm]' );
h = plot_spectrum( chemshift,spect_substrate ); title( label_disp ); hold on; xlabel( 'chemical shift [ppm]' );
hold on


    substrate_scale = 1;
map_substrate = squeeze( sum( abs( data( idx_substrate,:,: )),1 ));
