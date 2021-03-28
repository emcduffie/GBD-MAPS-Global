%FUNCTION: interpfn 
% interpolation of data (in) at any global grid (glat, glon) to  another
% global grid (lats, lons)
%11/20/2018
%--------------------------
%INPUTS: 
%in - original data (needs to be in lat, lon)
%glat, glon - original lat and lon of data
%lats,lons - new lat and lon of data
%OUTPUTS
%out - data gridded to the new lat and lon
%--------------------------

function out = interpfn_xMAPS(glat, glon, in, lats, lons)

[GLAT, GLON] = meshgrid(glat,glon);
[LATS, LONS] = meshgrid(lats,lons);

out =  interp2(GLAT,GLON,in',LATS,LONS,'linear')';