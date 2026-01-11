
%% Calculate Mutual informaiton of two variables
smoothingValue = 0.5;
nRep = 100;
nObs = 10000;
cols = colororder;

figure;
for k = 1 : nRep
    x = randn(nObs,1);
    y = randn(nObs,1);
    MI(k) = getMI(x,y,smoothingValue);
end

histogram(MI,FaceColor=cols(1,:)); hold on;

smoothingValue = 0.5;
nRep = 100;

for k = 1 : nRep
    x = randn(nObs,1);
    y = 2*x;
    MI(k) = getMI(x,y,smoothingValue);
end

histogram(MI,FaceColor=cols(2,:))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MI = getMI(x,y,smoothingValue)
% Prepare for histogram
nBins = 10;
binEdgesX = linspace(min(x),max(x),nBins+1);
binEdgesY = linspace(min(y),max(y),nBins+1);

% Histograms of X
hX = histogram(x,nBins,"BinEdges",binEdgesX,"Visible","off"); hold on;

% Histogram of Y
hY = histogram(y,nBins,"BinEdges",binEdgesY,"Visible","off");

% Histogram of XY
hXY = histogram2(x,y,nBins,XBinLimits=[binEdgesX(1), binEdgesX(end)],...
    YBinLimits=[binEdgesY(1), binEdgesY(end)],Visible="off");

% Add counts to smooth the histogram
XValues = hX.Values' + smoothingValue;
YValues = hY.Values' + smoothingValue;
XYValues = hXY.Values' + smoothingValue;

% Calculate probabilities
pX = XValues / sum(XValues);
pY = YValues / sum(YValues);
pXY = XYValues / sum(sum(XYValues));

% Prepare Grids for MI calculation
XGrid = repelem(pX,1,nBins);
YGrid = repelem(pY',nBins,1);

% get MI
PInd = XGrid.*YGrid;
MI = sum(sum(pXY.*log(pXY./PInd)));
end