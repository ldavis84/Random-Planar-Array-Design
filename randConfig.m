function [r,az] = randConfig(rFix,azFix,nRand,xGrid,yGrid,iAllow,azAllow)
% [r,az] = randConfig(rFix,azFix,nRand,xGrid,yGrid,iAllow,azAllow);
% Returns an allowed planar array of apertures.
% rFix        -- 3 x nfix positions of fixed apertures
% azFix       -- 1 x nfix azimuth orientations of the fixed apertures, deg
% nRand       -- number of random apertures to place
% xGrid,yGrid -- nx x ny X and Y positions of the grid points. xGrid(i,:)
%                is constant, yGrid(:,i) is constant
% iAllow      -- nx x ny logical array of which grid points are allowed
% azAllow     -- nx x ny x 2 allowed interval of azimuths for each node
%                [azAllow(i,j,1) azAllow(i,j,2)] is the interval in deg
%                Not necessarily in any particular branch cut.  E.G.,
%                suppose 10-20 deg is not allowed.  Put
%                azAllow(i,j,1:2) = [20 10+360];  Default is [0 360];
% r           -- 3 x nfix+nrand coordinates of the array apertures
% az          -- 1 x nfix+nrand azimuth rotations of the array apertures
%                in deg, in (-180,180]
[nx,ny] = size(xGrid);

if (~exist('azAllow','var') || isempty(azAllow))
    azAllow = zeros(nx,ny,2);
    azAllow(:,:,2) = 360;
end

iAllow = iAllow(:);
xGrid = xGrid(iAllow);
yGrid = yGrid(iAllow);

azAllow = reshape(azAllow,nx*ny,2);
azAllow = azAllow(iAllow,:);

nAllow = length(xGrid);

% pick nRand out of nAllow at random.

chosen = pickiofn(nRand,nAllow);

azAllow = azAllow(chosen,:);

% get configuration

rrand = [xGrid(chosen), yGrid(chosen), zeros(nRand,1)]';
r = [rFix, rrand];

azrand = azAllow(:,1) + rand(size(azAllow,1),1) .* (azAllow*[-1; 1]);

az = [azFix, azrand'];
az = to180(az);

end