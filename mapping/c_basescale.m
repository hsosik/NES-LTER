function s=c_basescale(lat)
% C_BASESCALE Rescale Figure 
%	C_BASESCALE, without any arguments, rescales the
%       current figure (where x is lon and y is lat), so
%       that dx=dy. It uses the mean latitude of the figure.
%
%	C_BASESCALE(LAT) uses LAT instead
%
%	S=C_BASESCALE gives the ratio between dlon/dlat at that 
%	latitude, which is used to rescale the map
%
%       Note this only makes sense if the area being plotted is
% 	quite small. Otherwise you should use Rich Pawlowicz' m_map,
%	for example. This requires sw_dist (from the CSIRO SeaWater
%	routines).
%
%	@Carlos Moffat, 2004

% ChangeLog:
% v1.0 05/10/06 Cleaned the code a bit

if exist('sw_dist','file')~=2
error('Err:Err','You need sw_dist, which can be found in the CSIRO SeaWater ToolKit \n routines, at http://woodshole.er.usgs.gov/operations/sea-mat/')
end

if nargin==0
    lat=mean(ylim);
end

a = sw_dist([0 1],[0 0],'km');
b = sw_dist([lat lat],[0 1],'km');
s=b/a;
set(gca,'DataAspectratio',[1 s 1]);
