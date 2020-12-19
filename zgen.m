% 8/2012, function out = zgen(domain, order)
%
% This function (has been long coming) creates a zernike polynomial over
% the domain specified by DOMAIN at an arbitrary order ORDER.  If DOMAIN is
% a single number, it is taken to be a pinhole of resolution specified by
% that number.  Flag specifies:
%
%     pt: computes OUT to be a list of values specified by coordinate pairs
%     specified by a N x 2 vector inputed as DOMAIN
% 
%     fn: returns a handle to a function to compute zernikes over DOMAIN

function out = zgen(domain, order, flag)
if nargin == 2
    flag = 'none';
end

switch flag
    case 'pt'
        xind = domain(:,1);
        yind = domain(:,2);
        [TH, R] = cart2pol(xind, yind);
        [n, m] = j2nm(order);
        Z = ZgenNM(n, m);
        out = Z(R, TH);
        out(R>1) = nan;
        ZgenNMshow(n, m, R, TH);
    case 'fn'
        [n, m] = j2nm(order);
        zRTH = ZgenNM(n,m);
        out = @(x,y) zRTH(sqrt(x.^2 + y.^2), atan2(y,x));
    case 'fnr'
        [n, m] = j2nm(order);
        zRTH = ZgenNM(n,m);
        out = @(r,th) zRTH(r, th);
        
    case 'none'
        if length(domain) == 1
            % assume DOMAIN is specifying resolution
            domain = pinhole(domain(1));
        end
        [sy, sx] = size(domain);
        xind = linspace(-1, 1, sx);
        yind = linspace(-1, 1, sy);
        [X Y] = meshgrid(xind, yind);
        idx = find(domain);
        [TH R] = cart2pol(X(idx), Y(idx));
        [n, m] = j2nm(order);
        Z = ZgenNM(n, m);
        z = Z(R, TH);
        out = domain;
        out(idx) = z;
    case 'basis'
        [sy, sx] = size(domain);
        xind = linspace(-1, 1, sx);
        yind = linspace(-1, 1, sy);
        [X, Y] = meshgrid(xind, yind);
        out = zeros(sx*sy, order + 1);
        out(:,1) = domain(:);
        for k = 1:order
            [n, m] = j2nm(k);
            zRTH = ZgenNM(n,m);
            out(:, k+1) = zRTH(sqrt(X(:).^2 + Y(:).^2), atan2(Y(:),X(:))).*domain(:);
        end
end


function Z = ZgenNM(n,m)
R = RgenNM(n,abs(m));
A = AgenNM(m);
Z = @(r, phi) R(r).*A(phi);


function ZgenNMshow(n,m, r, phi)
R = RgenNM(n,abs(m));
A = AgenNM(m);


% Forms radial function based on n,m
function f = RgenNM(n,m)
p = (n-m)/2;
f = @(r) 0;

for k = 0:p
    %fprintf('[%d, %d], [%d, %d]\n',n-k, k, n-2*k,  (n-m)/2-k);
    f = @(r) f(r) +  (-1)^k*nchoosek(n-k,k)*nchoosek(n-2*k,(n-m)/2-k)*r.^(n - 2*k);
end

% Forms azimuthal function based on n,m
function g = AgenNM(m)
if m > 0
    g = @(phi) cos(m*phi);
elseif m < 0
    g = @(phi) -sin(m*phi);
else
    g = @(phi) 1;
end

function [n,m] = j2nm(j)
% j = j+1; % 0 is piston
smct = 2;
ct = 0;
numInSmct = 3;
while j > 0
    smct = smct + 2;
    j = j - numInSmct;
    numInSmct = smct + 1;
    ct = ct + 1;
end
n = 2*ct;
m = 0;

for k = 1:abs(j)
    if isodd(k)
        n = n - 1;
        m = m + 1;
    end
end
if isodd(abs(j))
    m = -m;
end