function canSee = seePastBlock(r,dir,rectBlock)
% canSee = seePastBlock(r,dir,rectBlock);
% Checks whether positions r(:,i) can see past rectangular blockages in
% direction dir.  Note that, if a point is lying on the blockage it will
% not be blocked if the surface normal has the same sign as the direction.
% r         -- 3 x nr point positions to check
% dir       -- 3 x 1 unit look direction vector from the point
% recBlock  -- 1 x m cell array of structures convenient for specifying
%              m rectangles.  Elements:
% .ro       -- 3 x 1 center of rectangle in global coords
% .Qrg      -- 3 x 3 direction cosine matrix converting global to rectangle
%              coordinates
% .wx, .wy  -- width of rectangle in x, y
% cansee    -- 1 x nr logical vector, true if point can see past all blocks

m = length(rectBlock);
nr = size(r,2);
disttol = 10*eps;

canSee = true(1,nr);

if isempty(rectBlock)
    return;
end

for i = 1:m
    
    Qrg = rectBlock{i}.Qrg;
    wx = rectBlock{i}.wx;
    wy = rectBlock{i}.wy;
    ro = rectBlock{i}.ro;
    dr = Qrg * dir;
    
    ro = ro * ones(1,nr);
    rr = Qrg*(r - ro);
    
    if (abs(dr(3)) >= disttol)
        t = -rr(3,:) / dr(3);
        
        rint = rr + dr * t;
        inplane = abs(t) <= disttol;
        
        thisBlocked = (inplane & dr(3) <= 0) ...
            | (~inplane & abs(rint(1,:)) <= wx/2 ...
            & abs(rint(2,:)) <= wy/2 & t > 0);
        canSee = canSee & ~thisBlocked;
        
    else  % look direction in plane of blockage--canSee if it is skew
        
        canSee = canSee & abs(dr'*rr) > disttol;
    end
end

end