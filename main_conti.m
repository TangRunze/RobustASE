
clear all
close all

%% --- Parameters Setting ---
epsilon = 0.2;
epsilonInB = 0.01;
scaleB = 1;
scaleVar = 10;
m = 10;

nVertex = 150;
nBlock = 3;
theta = ones(1, nBlock);
dimLatentPosition = nBlock;
rho = repmat(1/nBlock, 1, nBlock);
iStart = 1;
iEnd = 100;



% block probability matrix
B = (0.5 - epsilonInB)*ones(nBlock, nBlock) + 2*epsilonInB*eye(nBlock);
B = scaleB^2*B;

%% --- Parallel Computing ---
% delete(gcp('nocreate'))
% parpool(nCore);


%% Monte Carlo Simulation

errorRateMean = zeros(1, iEnd);
errorRateHL = zeros(1, iEnd);
ADiffMean = zeros(1, iEnd);
ADiffHL = zeros(1, iEnd);
AXDiffMean = zeros(1, iEnd);
AXDiffHL = zeros(1, iEnd);

for iIter = iStart:iEnd
    adjMatrixTotal = zeros(m, nVertex, nVertex);
    adjMatrixSum = zeros(nVertex, nVertex);
    
    % Generate block assignment.
    tauStarInd = mnrnd(1, rho, nVertex);
    tauStar = [1:nBlock]*tauStarInd';
    
    for iGraph = ((iIter - 1)*m+1):(iIter*m)
        % Generate data if there does not exist one, otherwise read the
        % existing data.
        [adjMatrix, ~] = datagenerator_conti(nVertex, ...
            B, tauStar, epsilon, scaleVar, iGraph);
        adjMatrixTotal(iGraph - (iIter - 1)*m, :, :) = adjMatrix;
        adjMatrixSum = adjMatrixSum + adjMatrix;
    end
    
    adjMatrixMean = adjMatrixSum/m;
    adjMatrixHL = hlcalculator(adjMatrixTotal, m);
    
    PStar = B(tauStar, tauStar);
    
    % mean
    ADiffMean(iIter) = norm(adjMatrixMean - PStar);
    
    xHatMean = asge(adjMatrixMean, dimLatentPosition);
    gm = fitgmdist(xHatMean, nBlock, 'Replicates', 10);
    tauHatMean = cluster(gm, xHatMean)';
    errorRateMean(iIter) = errorratecalculator(tauStar, tauHatMean, nBlock);
    
    AXDiffMean(iIter) = norm(xHatMean*xHatMean' - PStar);
    
    % HL
    ADiffHL(iIter) = norm(adjMatrixHL - PStar);
    
    xHatHL = asge(adjMatrixHL, dimLatentPosition);
    gm = fitgmdist(xHatHL, nBlock, 'Replicates', 10);
    tauHatHL = cluster(gm, xHatHL)';
    errorRateHL(iIter) = errorratecalculator(tauStar, tauHatHL, nBlock);
    
    AXDiffHL(iIter) = norm(xHatHL*xHatHL' - PStar);
end

[mean(errorRateMean), mean(errorRateHL)]

% 1-sided sign-test HA errorRateMean > errorRateHL
tmpStats = sum(errorRateMean > errorRateHL);
pValueError = 1 - binocdf(tmpStats - 1, iEnd, 0.5)

[mean(ADiffMean), mean(ADiffHL)]

% 1-sided sign-test
tmpStats = sum(ADiffMean > ADiffHL);
pValueADiff = 1 - binocdf(tmpStats - 1, iEnd, 0.5)

[mean(AXDiffMean), mean(AXDiffHL)]

% 1-sided sign-test
tmpStats = sum(AXDiffMean > AXDiffHL);
pValueAXDiff = 1 - binocdf(tmpStats - 1, iEnd, 0.5)




%% --- Plot ---
% plot3(xHatMean(:, 1), xHatMean(:, 2), xHatMean(:, 3), '.')
% hold on;
% plot3(xHatHL(:, 1), xHatHL(:, 2), xHatHL(:, 3), 'r.')

%% --- Close Parallel Computing ---
% delete(gcp('nocreate'))
