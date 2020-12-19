
% This script adds programmatic zernikes to the 2-ray method to form a
% sheared basis and compares the condition numbers of various bases.
clc
clear all
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
Nz = 37; % number of zernikes
mask=logical(pinhole(N));
Ne=sum(mask(:));
NI = zeros(2*Ne,1); % null interferogram
betaBasis = zeros(2*Ne, Nz);
for k = 0:Nz

    if k == 0
        zernCouples = [0, 0];
 
    else
        zernCouples = [k, 1];
    end

    fprintf('Computing basis %d of %d\n', k, Nz);
    % Computes the opd from +1 and 0 grating orders
    [opds2x, opds2y] = twoRaySystem(X, Y, -Z,z1_um, z2s_cm*1e4, lambda_um, T_um, [0, 1], ...
            alpha, beta, zernCouples', NA);
    [opds2sx, opds2sy] = twoRaySystem(X, Y, -Z, z1_um, z2s_cm*1e4, lambda_um, T_um, [0, -1], ...
            alpha, beta, zernCouples', NA);
        
    opds2x=(-opds2x+ opds2sx)/2 .* pinhole(N);
    opds2y=(-opds2y+ opds2sy)/2 .* pinhole(N);
   % figure(1),mesh(opds2x)
        

    % Generate bases
    if k == 0
        NI = [opds2x(mask) ; opds2y(mask)]/lambda_um;  % Null interferogram
    else
        betaBasis(:,k) =  [opds2x(mask) ; opds2y(mask)]/lambda_um - NI;
    end
        
end

% Generate alternative bases
betaBasisNI = [betaBasis, NI];
gammaBasis = [NI, betaBasis + NI*ones(1,37)];

% Normalize bases:
betaBasisNorm = betaBasis./   sqrt((ones(2*Ne, 1)*sum(betaBasis.^2)));
betaBasisNINorm = betaBasisNI./   sqrt((ones(2*Ne, 1)*sum(betaBasisNI.^2)));
gammaBasisNorm = gammaBasis./   sqrt((ones(2*Ne, 1)*sum(gammaBasis.^2)));
%%
% Compute condition numbers
fprintf('\nCondition number of the beta basis is %0.3f\n', ...
            max(svd(betaBasisNorm))/min(svd(betaBasisNorm)));
fprintf('Condition number of the beta basis with NI as a vector is %0.3f\n', ...
            max(svd(betaBasisNINorm))/min(svd(betaBasisNINorm)));
fprintf('Condition number of the gamma basis is %0.3f\n', ...
            max(svd(gammaBasisNorm))/min(svd(gammaBasisNorm)));

% Decompose NI as a set of beta basis vectors:
coef = pinv(betaBasis)*NI;
residual = (sum(ones(2*Ne,1)*coef'.*betaBasis, 2) - NI);
fprintf('\nRMS error of representing NI as a linear combination of %d sheared zernikes: %0.5g waves\n', Nz, std(residual(:)));

sigZTerms = find(abs(coef) > 0.001);
fprintf('\nZernike terms that are larger than 1 mWave:\n------------------------------------------\n')
for k = 1:length(sigZTerms)
    fprintf('N = %2d: %6.1f mWaves\n', sigZTerms(k), 1000*coef(sigZTerms(k)));
end


%% 
% imagesc(abdom(pinhole(200), [0, 0, 0, 0, coef(4:end)']))
% title('Systematic aberration');