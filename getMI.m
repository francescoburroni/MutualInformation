function MI = getMI(x, y, nBins, smoothingValue, units, options)
% getMI - Calculate mutual information between two variables
%
% Syntax: MI = getMI(x, y, smoothingValue, units)
%
% Inputs:
%   x              - Numeric vector of observations for variable X
%   y              - Numeric vector of observations for variable Y
%   nBins          - Positive scalar for bins used in histogram calculations
%   smoothingValue - Positive scalar for Laplace smoothing to avoid log(0)
%                    (default: 0.5)
%   units          - String specifying units: 'nats' or 'bits'
%                    (default: 'nats')
%
% Output:
%   MI - Mutual information in specified units
%
% Description:
%   Computes the mutual information MI(X;Y) using histogram-based estimation
%   with 10 bins. Smoothing is applied to avoid zero probabilities.
%   Units can be 'nats' (natural log) or 'bits' (log base 2).
%
% Examples:
%   x = randn(100,1);
%   y = x.^2;
%   MI_nats = getMI(x, y, 0.5, 'nats');
%   MI_bits = getMI(x, y, 0.5, 'bits');

% Input validation
arguments
    x (:,1) double {mustBeNumeric, mustBeFinite, mustBeNonempty}
    y (:,1) double {mustBeNumeric, mustBeFinite, mustBeNonempty}
    nBins (1,1) double {mustBeInteger}
    smoothingValue (1,1) double {mustBeNumeric, mustBePositive, mustBeFinite} = 0.5
    units (1,1) string {mustBeMember(units, ["nats", "bits"])} = "bits"
    options.doPlot logical = false
    options.colorMap (1,1) string = "gray"
end

% Check that x and y have the same length
if length(x) ~= length(y)
    error('getMI:DimensionMismatch', 'Input vectors x and y must have the same length.');
end

% Check minimum bin number
if nBins < 2
    error('getMI:nBinsTooSmall', 'The number of bins (n=%d) has to be greater than 2.', nBins);
end

% Check minimum sample size
if length(x) < 10
    warning('getMI:SmallSampleSize', 'Sample size is very small (n=%d). Results may be unreliable.', length(x));
end

% Create bin edges for X and Y spanning their respective ranges
binEdgesX = linspace(min(x), max(x), nBins + 1);
binEdgesY = linspace(min(y), max(y), nBins + 1);

% Compute joint histogram for (X,Y)
hXY = histogram2(x, y, nBins, ...
    XBinLimits=[binEdgesX(1), binEdgesX(end)], ...
    YBinLimits=[binEdgesY(1), binEdgesY(end)], ...
    Visible="off",FaceColor="flat");

% Apply Laplace smoothing to histogram counts to avoid zero probabilities
XY = hXY.Values + smoothingValue;
pXY = XY / sum(XY, "all");

% Derive marginals from pXY
pX = sum(pXY, 1)';  % P(X) - marginal distribution
pY = sum(pXY, 2)';  % P(Y) - marginal distribution

% Create grids for computing independence assumption P(X)P(Y)
% XGrid(i,j) = P(X=i) for all j
XGrid = repelem(pX, 1, nBins);
% YGrid(i,j) = P(Y=j) for all i
YGrid = repelem(pY, nBins, 1);

% Compute product of marginals (independence assumption)
PInd = XGrid .* YGrid;  % PInd(i,j) = P(X=i) * P(Y=j)

% Select logarithm function based on desired units
if units == "bits"
    logFunc = @log2;  % Base-2 logarithm for bits
else
    logFunc = @log;   % Natural logarithm for nats
end

% Calculate mutual information using the formula:
% MI(X;Y) = sum_i sum_j P(X=i,Y=j) * log(P(X=i,Y=j) / (P(X=i)*P(Y=j)))
MI = sum(sum(pXY .* logFunc(pXY ./ PInd)));

% Plot

if options.doPlot 
   
    hXY.Normalization = "probability";
    hXY.Visible = "on";    
    colorbar
    colormap(options.colorMap)

end
end