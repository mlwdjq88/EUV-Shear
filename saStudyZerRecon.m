%%%
% This script adds programmatic zernikes to the 2-ray method to form a
% sheared basis and compares the condition numbers of various bases.

simMode = 'MET5';

switch simMode
    case 'MET5'
        lambda_um   = 0.0135;
        NA          = 0.5;
        T_um        = .234;
        z1_um       = 1;
        z2s_cm      = 2.2+1*z1_um/1e4; % Realistic numbers for detector

    case 'SERM'
        lambda_um   = 0.01356;
        NA          = 0.0875;
        T_um        = 1.544;
        z1_um       = 352;
        z2s_cm      = 12.7+1*z1_um/1e4; % Realistic numbers for detector
end

N           = 128;
z2_um       = z2s_cm * 1e4-z1_um;
halfD       = tan(asin(NA))*(z2_um);
xIdx        = linspace(-halfD, halfD, N);
yIdx        = linspace(-halfD, halfD, N);
[X, Y]      = meshgrid(xIdx,yIdx);

% Tilt angles
alpha   = 0;
beta    = 0;
eta     = 0;
gamma   = 0;

Z = X/cos(gamma)*sin(gamma)+Y/cos(eta)*sin(eta);

% set aberrations:
Nz = 1; 
mask=logical(pinhole(N));
Ne=sum(mask(:));
Aberration = zeros(2*Ne, Nz);
zernCouples = [1, 0.1;2, 0.2;3 0.4;4, 0.3;5, 0.2;6, 0.4;7, 0.2;8, 0.3;9, 0.1];

[opds2x, opds2y] = twoRaySystem(X, Y, -Z,z1_um, z2s_cm*1e4, lambda_um, T_um, [0, 1], ...
        alpha, beta, zernCouples', NA);
[opds2sx, opds2sy] = twoRaySystem(X, Y, -Z, z1_um, z2s_cm*1e4, lambda_um, T_um, [0, -1], ...
        alpha, beta, zernCouples', NA);

opds2x=(-opds2x+ opds2sx)/2 .* pinhole(N);
opds2y=(-opds2y+ opds2sy)/2 .* pinhole(N);
Aberration(:,1) =  [opds2x(mask) ; opds2y(mask)]/lambda_um - NI;
coefA = pinv(betaBasis)*Aberration;
