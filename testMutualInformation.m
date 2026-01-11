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
