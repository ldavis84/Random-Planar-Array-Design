% test seePastBlock.m
clear;

h = 1;   % height of the wall
w = 2;   % half-width

% four walls and a floor

abcd(:,:,1) = repmat([w;w;h],1,4) .* [1 1 1 1; 1 1 -1 -1; 0 1 1 0];
abcd(:,:,2) = repmat([w;w;h],1,4) .* [-1 -1 -1 -1; 1 1 -1 -1; 0 1 1 0];
abcd(:,:,3) = repmat([w;w;h],1,4) .* [1 1 -1 -1; 1 1 1 1; 0 1 1 0];
abcd(:,:,4) = repmat([w;w;h],1,4) .* [1 1 -1 -1; -1 -1 -1 -1; 0 1 1 0];
abcd(:,:,5) = repmat([w;w;h],1,4) .* [1 -1 -1 1; 1 1 -1 -1; 0 0 0 0];

rec = rec2struct(abcd);

% vantage point

r = [0;0;0];

az = 135;
el = 19.472;
dir = [cosd(az).*cosd(el); sind(az).*cosd(el); sind(el)];

canSee = seePastBlock(r,dir,rec)