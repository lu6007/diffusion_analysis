function [f,noise] = region_wiener2(g, nhood, BW, type)
%WIENER2 Perform 2-D adaptive noise-removal filtering.
%   WIENER2 lowpass filters an intensity image that has been degraded by
%   constant power additive noise. WIENER2 uses a pixel-wise adaptive Wiener
%   method based on statistics estimated from a local neighborhood of each
%   pixel.
%
%   [J,NOISE] = WIENER2(I,[M N]) also estimates the additive noise power
%   before doing the filtering. WIENER2 returns this estimate as NOISE.
%
%   Class Support
%   -------------
%   The input image I can be uint8, uint16, int16, double, or single.  The
%   output image J has the same class as I.
%
%   Example
%   -------
%       RGB = imread('saturn.png');
%       I = rgb2gray(RGB);
%       J = imnoise(I,'gaussian',0,0.005);
%       K = wiener2(J,[5 5]);
%       figure, imshow(J), figure, imshow(K)
%
%   See also FILTER2, MEDFILT2.

%   modified based on the adaptive linear filer winer2()

% Reference: "Two-Dimensional Signal and Image Processing" by 
% Jae S. Lim, p. 548, equations 9.44 - 9.46.

classin = class(g);
classChanged = false;
if ~isa(g, 'double')
  classChanged = true;
  g = im2double(g);
end



% 
region = filter2(ones(nhood),BW);
region = region + ~BW; % avoid divide by zero

%
g_region = g.*double(BW);

% Estimate the local mean of f.
localMean_region = filter2(ones(nhood), g_region) ./region;

% Estimate of the local variance of f.
localVar_region = filter2(ones(nhood), g_region.^2) ./ region ...
    - localMean_region.^2;


% Estimate the noise power if necessary.
localMean = filter2(ones(nhood), g) / prod(nhood);
localVar = filter2(ones(nhood), g.^2) / prod(nhood) - localMean.^2;
noise = mean2(localVar);

% Compute result
% f = localMean + (max(0, localVar - noise) ./ ...
%           max(localVar, noise)) .* (g - localMean);
%
% Computation is split up to minimize use of memory
% for temp arrays.
f = g_region - localMean_region;
g_region = localVar_region - noise; 
g_region = max(g_region, 0);
localVar_region = max(localVar_region, noise);
localVar_region = localVar_region+~BW; % avoid divinding by zero
f = f ./ localVar_region;
f = f .* g_region;
f = f + localMean_region;
f = f.*BW;

% if not unint16, us changeclass(), see wiener2()
if classChanged
  f = im2uint16(f);
end





