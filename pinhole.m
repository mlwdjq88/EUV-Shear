function p = pinhole(D, r, c)

if nargin == 1
    r = D;
    c = D;
end

if D == 1
    p = pad2(1, r,c);
    return;
elseif D == 2
    p = pad2(1, r,c);
    return;
end

if nargin == 2
    p = pinhole(min(D,r), D, r);
    return;
end

[X,Y] = meshgrid(  linspace(-1,1,D)  );
frame = zeros(D);
frame(X.^2+Y.^2 <= 1) = 1;

p = pad2(frame, r,c);