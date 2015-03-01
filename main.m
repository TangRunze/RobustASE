
%% --- Parameters Setting ---
epsilon = 0.1;
m = 10;
nVertex = 150;
nBlock = 3;
epsilonInB = 0.1;
theta = ones(1, nBlock);
dimLatentPosition = nBlock;
rho = repmat(1/nBlock, 1, nBlock);
iStart = 1;
iEnd = 10;

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

for iIter = iStart:iEnd
    adjMatrixTmp = zeros(nVertex, nVertex);
    for iGraph = ((iIter-1)*m+1):iIter*m
        % Generate data if there does not exist one, otherwise read the
        % existing data.
        [adjMatrix, tauStar, xStar] = datagenerator(nVertex, nBlock, ...
            dimLatentPosition, B, nuStar, rho, epsilon, iGraph);
        adjMatrixTmp = adjMatrixTmp + adjMatrix;
    end
    adjMatrixMean = adjMatrixTmp/m;
    adjMatrixHL = hlcalculator(adjMatrixTmp, m);
    
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




%% --- Close Parallel Computing ---
% delete(gcp('nocreate'))
