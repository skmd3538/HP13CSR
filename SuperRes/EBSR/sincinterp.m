function Zrecon = sincinterp ( Zsampl , Width , Height)
% Zsampl i s the sampled s i g n a l
% Width  i s  the  d e s i r e d width  of  the  recon s tru cted	s i g n a l
% Height  i s  the  d e s i r e d height  of  the  recon s tru cted	s i g n a l
%  Zrecon  i s  the recon s tru cted	s i g n a l
[ Nx, Ny, Nz] = size ( Zsampl) ;

Zrecon = zeros ( Height , Width, Nz) ;
x=linspace (-Nx/2 ,Nx/2 , Width ) ;
y=linspace (-Ny/2 ,Ny/2 , Height ) ;
for d = 1: Nz
    Zr = zeros( Height , Width);
    for n  =  1 : Nx
        for m =  1 : Ny
            den = pi * ( x - ( n-0.5*Nx ) ) ;
            num = sin ( den ) ;
            ind = find ( den == 0 ) ;
            if ~isempty ( ind )
                den ( ind ) = 1; num( ind ) = 1;
            end
            sincx = num./den ;
            sincx  =  repmat ( sincx ,	[ Height ,	1 ] ) ;
            den = pi * ( y - (m-0.5*Ny ) ) ;
            num = sin ( den ) ;
            ind = find ( den == 0 ) ;
            if ~isempty ( ind )
                den ( ind ) = 1 ; num( ind ) = 1 ;
            end
            sincy = num./den ;
            sincy = repmat ( sincy' , [ 1 , Width ] ) ;
            Zr = Zr + double(Zsampl(m, n, d)) .* sincx .* sincy ;
        end
    end
    Zrecon(:,:,d) = Zr;
end
end
