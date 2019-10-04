function Main()

% setenv('PATH', [getenv('PATH') ':/Library/TeX/texbin/:/usr/local/bin/']);
% close all; clear all; clc
% format LONG

%% READ THE ORIGINAL IMAGE
Im          = 255*imread('OD_05_ThCleaned.png');
Im          = double(imresize(Im, 0.5)); 

%% IMAGE PROCESSING STEP

% Transform the image to grayscale
[h, w, c]   = size(Im);
if (c>1), I0 = rgb2gray(Im); else, I0 = Im; end

imwrite(I0, 'CContour11.png', 'png')

% Normalize to [0,1]
In          = Normalize(double(I0), 1, 0);

% Median filtering 5 x 5 or 9 x 9
Imf         = medfilt2(In, [5 5],'symmetric'); 

% Low-pass filtering via Gaussian with mu = 0, sigma = 1, win_size = 5;
If          = I2Grad(Imf, 0, 1, 5);

% Image binarization: points extraction and eliminate border artifacts
Dborder     = 5;
Ix = If(Dborder : h - Dborder, Dborder : w - Dborder);
Ix(1,1)     = 1; 
Ix(size(Ix))= 1; 

% Show current state
figure('Name','Binarized Image'), imshow(Ix, []), 

%% DATA POINTS EXTRACTION

% Find all data points
[Py, Px]    = find(Ix>=1);
Py          = [Py ; min(Py) ; max(Py)];  
Py          = h - Py + 1;

Px          = [Px ; min(Px) ; max(Py)]; 
%Px          = w - Px + 1;

% Show current state
figure('Name','Extracted Data Points'), 
plot(Px, Py, '.b', 'MarkerSize', 2, 'LineWidth', 2)

%% CONIC OPTIMIZATION

% Initialize variables
Results     = [];
BestR       = Inf;
MPb         = 1;
BParab      = [];
th          = 6; %2,  0.005*min(h,w) % Inliers ~ 1 %  th < 10
Gk          = 0;

% Find the optical conic
[Parabola1,  niterA] = EvolDiff(Px, Py, th, Gk);

% Determine model Errors
[MAE, RMSE, MP, X1, Y1, X2, Y2] = ParabolaErrors(Parabola1, Px, Py, th);

% Show the achieved parabola
figure('Name','Achieved Parabola'), 
plotParabolaXY(Parabola1, [Px, Py], 3); pause(0.1)

% Show full result
figure('Name','Result');
[X1,Y1,X2,Y2] = plotParabolaXY(Parabola1, [Px, Py], 0); 
imshow(Im), axis on, hold on,
X1          = X1 + Dborder;         X2          = X2 + Dborder; 
Y1          = Y1 - Dborder;         Y2          = Y2 - Dborder;
plot (X1,h - Y1 + 1,'.b', X2, h - Y2 + 1, '.b', 'MarkerSize', 2, 'LineWidth', 2), 
pause(0.1), hold off
end

%% Normalize an image Im to [minValue,maxValue]
function ImOut = Normalize(Im, minValue, maxValue)
maxImValue  = max(Im(:)); 
minImValue  = min(Im(:));

ImOut       = (minValue - maxValue) / (maxImValue - minImValue) * ...
              (Im - minImValue) + maxValue;
end

%% Low-pass filtering via Gaussian with mu = 0, sigma = 1, win_size = 5;
function IoG = I2Grad(Im, mu, sigma, samples)
% Default values
if nargin < 2, mu = 0; sigma = 1; samples = 5; end

% Define x series
lowerLimit  = -3*sigma; 
upperLimit  = +3*sigma;
x           = linspace(lowerLimit, upperLimit, samples); 

% Gaussian function and its derivative
Gauss       = @(x, mu, sig) exp(-(x-mu).^2/2/sig^2)/sqrt(2*pi)/sig;
derGauss    = @(x, mu, sig) -(x-mu).*Gauss(x, mu, sig)/sqrt(2*pi)/sig^3;

% Find Gaussian samples and normalize them
DG          = derGauss(x,mu,sigma);
normDG      = DG / max(DG); 

% Filter the image
IoG         = sqrt(conv2(Im, normDG,'same').^2 + conv2(Im, normDG','same').^2);
end
