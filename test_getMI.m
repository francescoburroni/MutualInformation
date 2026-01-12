%% Calculate Mutual informaiton of two variables
smoothingValue = 0.5;
nBins = 25;
nRep = 100;
nObs = 10000;
cols = colororder;

figure;
for k = 1 : nRep
    x = randn(nObs,1);
    y = randn(nObs,1);
    MI(k) = getMI(x,y,nBins,smoothingValue);
end

histogram(MI,FaceColor=cols(1,:)); hold on;

smoothingValue = 0.5;
nRep = 100;

for k = 1 : nRep
    x = randn(nObs,1);
    y = 2*x;
    MI(k) = getMI(x,y,nBins,smoothingValue);
end

histogram(MI,FaceColor=cols(2,:))

%%

load("x1.csv")
load("y1.csv")
load("x2.csv")
load("y2.csv")

figure(Theme="light")
nBins = 20;
smoothingValue = 0.5;
MI1 = getMI(x1,y1,nBins,smoothingValue,"bits",doPlot=true,colorMap="bone");
fprintf("MI of x1 and y1 is %.3f\n",MI1)
MI2 = getMI(x2,y2,nBins,smoothingValue,"bits");
fprintf("MI of x2 and y2 is %.3f\n",MI2)