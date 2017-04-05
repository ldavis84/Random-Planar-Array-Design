%--------------------------------------------------------------------------
% Test optRandConfig tool.
%--------------------------------------------------------------------------

clear;
rng(2);

nTrials = 500;
lamEval = [10 100];

azTrue = [0 45 90 135 0 45 90 135];
elTrue = [30 30 30 30 80 80 80 80];

truePol = [0 0];

SNRdB = 0: 0.5 : 30;

defdB = 20;
verbose = true;

%--------------------------------------------------------------------------
% Azimuth and Elevation Points
%--------------------------------------------------------------------------

azvals = 0:1:360;
naz = length(azvals);
elvals = 0:1:88;
nel = length(elvals);

azGrid = repmat(azvals,nel,1);
elGrid = repmat(elvals',1,naz);

%--------------------------------------------------------------------------
% The rectangles blocking.
%--------------------------------------------------------------------------

hWall = 0.1;   % height of the wall
wWall = 8;   % half-width

% four walls and a floor

% abcd(:,:,1) = repmat([wWall;wWall;hWall],1,4)...
%     .* [1 1 1 1; 1 1 -1 -1; 0 1 1 0];
% abcd(:,:,2) = repmat([wWall;wWall;hWall],1,4)...
%     .* [-1 -1 -1 -1; 1 1 -1 -1; 0 1 1 0];
% abcd(:,:,3) = repmat([wWall;wWall;hWall],1,4)...
%     .* [1 1 -1 -1; 1 1 1 1; 0 1 1 0];
% abcd(:,:,4) = repmat([wWall;wWall;hWall],1,4)...
%     .* [1 1 -1 -1; -1 -1 -1 -1; 0 1 1 0];

abcd = [];

%--------------------------------------------------------------------------
% Define allowable grid and fixed apertures.
%--------------------------------------------------------------------------

dGrid = min(lamEval) / 5;
LGrid = 30;

xGrid = -LGrid/2 : dGrid : LGrid/2;
nGrid = length(xGrid);

xGrid = repmat(xGrid(:),1,nGrid);
yGrid = xGrid';

iAllow = sqrt(xGrid.^2 + yGrid.^2) >= 10;  % only stuff in the circle

rFix = LGrid/2* [-1 1 1 -1; 1 1 -1 -1; 0 0 0 0];   % four corners
azFix = 45 + 90*(0:3);

azAllow = [];   % no restrictions on placement of random stuff

nRand = 32 - length(azFix);    % how many random apertures allowed

%--------------------------------------------------------------------------
% Make the call.
%--------------------------------------------------------------------------

[rBest,azBest,azThBest,elThBest,hist] = ...
    optRandConfig(@horizB,nTrials,lamEval,azTrue,elTrue,abcd,...
    defdB,azGrid,elGrid,truePol,SNRdB,...
    rFix,azFix,nRand,xGrid,yGrid,iAllow,azAllow,verbose);

euler = [azBest; zeros(2,length(azBest))];

M = 5;
arrayAthley(...
    @horizB,SNRdB,lamEval(1),rBest,euler,M,[azTrue(2) elTrue(2)],...
    truePol,azGrid,elGrid,'ExampleBestRandom');

save(mfilename);