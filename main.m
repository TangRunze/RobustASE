
%% --- Parameters Setting ---
epsilon = 0.3;
m = 10;
nVertex = 150;
nBlock = 3;
epsilonInB = 0.07;
theta = ones(1, nBlock);
dimLatentPosition = nBlock;
rho = repmat(1/nBlock, 1, nBlock);
iStart = 1;
iEnd = 100;

% block probability matrix
B = (0.5 - epsilonInB)*ones(nBlock, nBlock) + 2*epsilonInB*eye(nBlock);

% true nu (nBlock-by-dimLatentPosition)
% nuStar = chol(B)';
[U, S, ~] = svds(B, dimLatentPosition);
nuStar = U*sqrt(S);

% true spectral graph embedding Xhat (nBlock-by-dimLatentPosition)
xHatHL = asge(B, dimLatentPosition);

%% --- Parallel Computing ---
% delete(gcp('nocreate'))
% parpool(nCore);


%% Monte Carlo Simulation

errorRateMean = zeros(1, iEnd);
errorRateHL = zeros(1, iEnd);
aMean = zeros(1, iEnd);
aHL = zeros(1, iEnd);

for iIter = iStart:iEnd
    
    % Generate data if there does not exist one, otherwise read the
    % existing data.
    [adjMatrixSum, tauStar] = datagenerator(nVertex, nBlock, ...
        dimLatentPosition, B, nuStar, rho, epsilon, m, iIter);
    
    adjMatrixMean = adjMatrixSum/m;
    adjMatrixHL = hlcalculator(adjMatrixSum, m);
    
    % mean
    xHatMean = asge(adjMatrixMean, dimLatentPosition);
    gm = fitgmdist(xHatMean, nBlock, 'Replicates', 10);
    tauHatMean = cluster(gm, xHatMean)';
    errorRateMean(iIter) = errorratecalculator(tauStar, tauHatMean, nBlock);
    
    % HL
    xHatHL = asge(adjMatrixHL, dimLatentPosition);
    gm = fitgmdist(xHatHL, nBlock, 'Replicates', 10);
    tauHatHL = cluster(gm, xHatHL)';
    errorRateHL(iIter) = errorratecalculator(tauStar, tauHatHL, nBlock);
    
end

errorRateMean

errorRateHL


%% --- Close Parallel Computing ---
% delete(gcp('nocreate'))
