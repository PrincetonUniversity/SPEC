function write_spec_Bxyz( data, filename )

nt = length(data.theta);
nz = length(data.phi);

fid = fopen( filename, 'w' );

fprintf( fid, 'theta           phi             x                y                z                Bx               By               Bz\n')


for it = 1:nt
  for iz = 1:nz
  
    theta = data.theta(it);
    phi   = data.phi(iz);
    x     = data.x(1,it,iz);
    y     = data.y(1,it,iz);
    z     = data.z(1,it,iz);
    Bx    = data.Bx(1,it,iz);
    By    = data.By(1,it,iz);
    Bz    = data.Bz(1,it,iz);

    
      %fprintf(fid,'%s', A{i});
    fprintf( fid, '%12.8E  %12.8E  %+12.8E  %+12.8E  %+12.8E  %+12.8E  %+12.8E  %+12.8E\n', ...
             theta, phi, x, y, z, Bx, By, Bz)

  end
end

fclose(fid);
end
