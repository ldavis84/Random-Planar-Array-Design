function angles = to180(anglesin)
% angles = to180(anglesin)
% This ensures that the input matrix of angles is mapped to the
% branch cut (-180,180].

[n,m] = size(anglesin);
onemat = ones(n,m);

kmat = ceil ( (anglesin - 180*onemat) / 360 );
angles = anglesin - kmat * 360;
