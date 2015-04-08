
clear all
close all

%% --- Parameters Setting ---

seed = 12345;
rng(seed);

epsilon = 0.2;
epsilonInB = 0.3;
scaleB = 1;
c = 20;
m = 10;
q = 0.90;

nVertex = 100;
nBlock = 2;
theta = ones(1, nBlock);
dimLatentPosition = nBlock;
rho = repmat(1/nBlock, 1, nBlock);
iStart = 1;
iEnd = 1;



% block probability matrix
B = (0.5 - epsilonInB)*ones(nBlock, nBlock) + 2*epsilonInB*eye(nBlock);
B = scaleB^2*B;

%% --- Parallel Computing ---
% delete(gcp('nocreate'))
% parpool(nCore);


%% Monte Carlo Simulation

options = optimoptions('fmincon', 'TolX', 1e-6, 'MaxIter', 10000, ...
    'MaxFunEvals', 10000);

errorRateMean = zeros(1, iEnd);
errorRateLq = zeros(1, iEnd);
ADiffMean = zeros(1, iEnd);
ADiffLq = zeros(1, iEnd);
AXDiffMean = zeros(1, iEnd);
AXDiffLq = zeros(1, iEnd);

for iIter = iStart:iEnd
    adjMatrixTotal = zeros(m, nVertex, nVertex);
    adjMatrixSum = zeros(nVertex, nVertex);
    adjMatrixLq = zeros(nVertex, nVertex);
    
    % Generate block assignment.
    tauStarInd = mnrnd(1, rho, nVertex);
    tauStar = (1:nBlock)*tauStarInd';
    
    for iGraph = ((iIter - 1)*m+1):(iIter*m)
        % Generate data if there does not exist one, otherwise read the
        % existing data.
        [adjMatrix, ~] = datagenerator_lq(nVertex, ...
            B, tauStar, epsilon, c, iGraph);
        adjMatrixTotal(iGraph - (iIter - 1)*m, :, :) = adjMatrix;
        adjMatrixSum = adjMatrixSum + adjMatrix;
    end
    
    adjMatrixMean = adjMatrixSum/m;
    for i = 1:nVertex
        for j = (i + 1):nVertex
            [i j]
%             adjMatrixLq(i, j) = fsolve(@(theta) ...
%                 objectivefun(theta, q, squeeze(adjMatrixTotal(:,i,j))), ...
%                 adjMatrixMean(i, j));
            
            adjMatrixLq(i, j) = lqsolve(squeeze(adjMatrixTotal(:,i,j)), q);
            
            adjMatrixLq(j, i) = adjMatrixLq(i, j);
        end
    end
    
    PStar = B(tauStar, tauStar);
    
    % Diagonal Augmentation.
%     adjMatrixMeanDA = adjMatrixMean + adjMatrixMean';
%     adjMatrixMeanDA = 1*(adjMatrixMeanDA > 0);
%     adjMatrixMeanDA = adjMatrixMeanDA + diag(sum(adjMatrixMeanDA))/...
%         (size(adjMatrixMeanDA, 1) - 1);
%     
%     adjMatrixLqDA = adjMatrixLq + adjMatrixLq';
%     adjMatrixLqDA = 1*(adjMatrixLqDA > 0);
%     adjMatrixLqDA = adjMatrixLqDA + diag(sum(adjMatrixLqDA))/...
%         (size(adjMatrixLqDA, 1) - 1);
    
    adjMatrixMeanDA = adjMatrixMean;
    adjMatrixLqDA = adjMatrixLq;

    % ASGE
    xHatMean = asge(adjMatrixMeanDA, dimLatentPosition);
    gmMean = fitgmdist(xHatMean, nBlock, 'Replicates', 10);
    tauHatMean = cluster(gmMean, xHatMean)';
    pTauHatMean = posterior(gmMean, xHatMean)';
    muHatMean = gmMean.mu;
    sigmaHatMean = gmMean.Sigma;
    
    xHatLq = asge(adjMatrixLqDA, dimLatentPosition);
    gmLq = fitgmdist(xHatLq, nBlock, 'Replicates', 10);
    tauHatLq = cluster(gmLq, xHatLq)';
    pTauHatLq = posterior(gmLq, xHatLq)';
    muHatLq = gmLq.mu;
    sigmaHatLq = gmLq.Sigma;
    
    scatterplot(xHatMean, tauHatMean, muHatMean, sigmaHatMean, 'r', 'Mean');
    hold on;
    scatterplot(xHatLq, tauHatLq, muHatLq, sigmaHatLq, 'b', 'Lq');
    plotLegend = legend('Mean: embedded points', ...
        'Mean: mean of cluster', 'Mean: 95% ellipse', ...
        'Lq: embedded points', 'Lq: mean of cluster', 'Lq: 95% ellipse');
    set(plotLegend, 'FontSize', 14);
    hold off;
    
    plotTitle = title(strcat('n = ', num2str(nVertex), ', K = ', ...
        num2str(nBlock), ', epsilon = ', num2str(epsilon), ', diag = ', ...
        num2str(0.5 + epsilonInB), ', offdiag = ', ...
        num2str(0.5 - epsilonInB), ', q = ', num2str(q), ...
        ', c = ', num2str(c), ', m = ', num2str(m), ...
        ', seed = ', num2str(seed)));
    set(plotTitle, 'FontSize', 14);
    
    % mean
    ADiffMean(iIter) = norm(adjMatrixMean - PStar);
    
%     xHatMean = asge(adjMatrixMean, dimLatentPosition);
%     gm = fitgmdist(xHatMean, nBlock, 'Replicates', 10);
%     tauHatMean = cluster(gm, xHatMean)';
    errorRateMean(iIter) = errorratecalculator(tauStar, tauHatMean, nBlock);
    
    AXDiffMean(iIter) = norm(xHatMean*xHatMean' - PStar);
    
    % Lq
    ADiffLq(iIter) = norm(adjMatrixLq - PStar);
    
%     xHatHL = asge(adjMatrixLq, dimLatentPosition);
%     gm = fitgmdist(xHatHL, nBlock, 'Replicates', 10);
%     tauHatHL = cluster(gm, xHatHL)';
    errorRateLq(iIter) = errorratecalculator(tauStar, tauHatLq, nBlock);
    
    AXDiffLq(iIter) = norm(xHatLq*xHatLq' - PStar);
end

[mean(errorRateMean), mean(errorRateLq)]

% 1-sided sign-test HA errorRateMean > errorRateHL
tmpStats = sum(errorRateMean > errorRateLq);
pValueError = 1 - binocdf(tmpStats - 1, iEnd, 0.5)

[mean(ADiffMean), mean(ADiffLq)]

% 1-sided sign-test
tmpStats = sum(ADiffMean > ADiffLq);
pValueADiff = 1 - binocdf(tmpStats - 1, iEnd, 0.5)

[mean(AXDiffMean), mean(AXDiffLq)]

% 1-sided sign-test
tmpStats = sum(AXDiffMean > AXDiffLq);
pValueAXDiff = 1 - binocdf(tmpStats - 1, iEnd, 0.5)




%% --- Plot ---
% plot3(xHatMean(:, 1), xHatMean(:, 2), xHatMean(:, 3), '.')
% hold on;
% plot3(xHatHL(:, 1), xHatHL(:, 2), xHatHL(:, 3), 'r.')

%% --- Close Parallel Computing ---
% delete(gcp('nocreate'))
