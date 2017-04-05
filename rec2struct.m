function rec = rec2struct(abcd)
% rec = rec2struct(abcd);
% Converts rectangle corners into structures.  No check on geometry.  Note
% the outward facing normal (local z) uses the right hand rule on the
% traversal of the vertices provided.
% abcd      -- 3 x 4 x m corners of rectangles, in order around edges
% rec       -- 1 x m cell array of structures convenient for specifying
%              m rectangles.  Elements:
% .ro       -- 3 x 1 center of rectangle in global coords
% .Qrg      -- 3 x 3 direction cosine matrix converting global to rectangle
%              coordinates
% .wx, .wy  -- width of rectangle in its own local x, y

if (isempty(abcd))
    rec = [];
    return;
end

m = size(abcd,3);
rec = cell(1,m);

ok = false(1,m);

for i = 1:m
    ro = abcd(:,1,i);
    dx = abcd(:,2,i) - ro;
    dy = abcd(:,4,i) - ro;
    
    wx = norm(dx);
    wy = norm(dy);
    
    ok(i) = (abs(wx) >= 100*eps && abs(wy) >= 100*eps);
    
    if ok(i)
        
        x = dx / wx;
        y = dy / wy;
        z = cross(x,y);
        
        ro = mean(abcd(:,:,i),2);
    
        Qrg = [x y z]';
        
        rec{i}.wx = wx;
        rec{i}.wy = wy;
        rec{i}.ro = ro;
        rec{i}.Qrg = Qrg;
    end
    
end

end
