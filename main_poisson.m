
clear all
close all

%% --- Parameters Setting ---

seed = 12345;
rng(seed);

epsilon = 0.2;
epsilonInB = 0.3;
scaleB = 1;
c = 20;
m = 20;
q = 0.90;

nVertex = 100;
nBlock = 2;
dimLatentPosition = nBlock;
rho = repmat(1/nBlock, 1, nBlock);
iStart = 1;
iEnd = 1;

% block probability matrix
B = (0.5 - epsilonInB)*ones(nBlock, nBlock) + 2*epsilonInB*eye(nBlock);
B = scaleB^2*B;
[U, S, V] = svd(B);
muStar = U*sqrt(S);

% plot(muStar(:, 1), muStar(:, 2), '.b');
% hold on;
theta = -3*pi/4;
rotationMatrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];
muStar = muStar*rotationMatrix';
% plot(muStar(:, 1), muStar(:, 2), '.r');
% hold off;

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
    adjMatrixTotal0 = zeros(m, nVertex, nVertex);
    adjMatrixSum0 = zeros(nVertex, nVertex);
    adjMatrixLq0 = zeros(nVertex, nVertex);

    % Generate block assignment.
    tauStarInd = mnrnd(1, rho, nVertex);
    tauStar = (1:nBlock)*tauStarInd';
    
    for iGraph = ((iIter - 1)*m+1):(iIter*m)
        % Generate data if there does not exist one, otherwise read the
        % existing data.
        [adjMatrix, ~, adjMatrix0] = datagenerator_poisson(nVertex, B, ...
            tauStar, epsilon, c, iGraph);
        adjMatrixTotal(iGraph - (iIter - 1)*m, :, :) = adjMatrix;
        adjMatrixSum = adjMatrixSum + adjMatrix;
        adjMatrixTotal0(iGraph - (iIter - 1)*m, :, :) = adjMatrix0;
        adjMatrixSum0 = adjMatrixSum0 + adjMatrix0;
    end
    
    adjMatrixMean = adjMatrixSum/m;
    adjMatrixMean0 = adjMatrixSum0/m;
    for i = 1:nVertex
        for j = (i + 1):nVertex
            [i j]
%             adjMatrixLq(i, j) = fsolve(@(theta) ...
%                 objectivefun(theta, q, squeeze(adjMatrixTotal(:,i,j))), ...
%                 adjMatrixMean(i, j));
            adjMatrixLq(i, j) = lqsolvepoisson(squeeze(adjMatrixTotal(:,i,j)), q);
            adjMatrixLq(j, i) = adjMatrixLq(i, j);
            adjMatrixLq0(i, j) = lqsolvepoisson(squeeze(adjMatrixTotal0(:,i,j)), q);
            adjMatrixLq0(j, i) = adjMatrixLq0(i, j);
        end
    end
    
    PStar = B(tauStar, tauStar);
    
    % Diagonal Augmentation.
    adjMatrixMeanDA = adjMatrixMean;
    adjMatrixMeanDA = adjMatrixMeanDA + diag(sum(adjMatrixMeanDA, 2))/...
        (size(adjMatrixMeanDA, 2) - 1);
    
    adjMatrixLqDA = adjMatrixLq;
    adjMatrixLqDA = adjMatrixLqDA + diag(sum(adjMatrixLqDA, 2))/...
        (size(adjMatrixLqDA, 2) - 1);
    
    adjMatrixMeanDA0 = adjMatrixMean0;
    adjMatrixMeanDA0 = adjMatrixMeanDA0 + diag(sum(adjMatrixMeanDA0, 2))/...
        (size(adjMatrixMeanDA0, 2) - 1);
    
    adjMatrixLqDA0 = adjMatrixLq0;
    adjMatrixLqDA0 = adjMatrixLqDA0 + diag(sum(adjMatrixLqDA0, 2))/...
        (size(adjMatrixLqDA0, 2) - 1);

    
    
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
    
    
    xHatMean0 = asge(adjMatrixMeanDA0, dimLatentPosition);
    gmMean0 = fitgmdist(xHatMean0, nBlock, 'Replicates', 10);
    tauHatMean0 = cluster(gmMean0, xHatMean0)';
    pTauHatMean0 = posterior(gmMean0, xHatMean0)';
    muHatMean0 = gmMean0.mu;
    sigmaHatMean0 = gmMean0.Sigma;
    
    xHatLq0 = asge(adjMatrixLqDA0, dimLatentPosition);
    gmLq0 = fitgmdist(xHatLq0, nBlock, 'Replicates', 10);
    tauHatLq0 = cluster(gmLq0, xHatLq0)';
    pTauHatLq0 = posterior(gmLq0, xHatLq0)';
    muHatLq0 = gmLq0.mu;
    sigmaHatLq0 = gmLq0.Sigma;
    
    
    % Rotate xHatMean to match muStar
    wMean = procrustes(muHatMean, muStar);
    muHatMean = muHatMean*wMean;
    xHatMean = xHatMean*wMean;
    for i = 1:nBlock
        sigmaHatMean(:, :, i) = wMean'*squeeze(sigmaHatMean(:, :, i))*wMean;
    end
    
    % Rotate xHatLq to match muStar
    wLq = procrustes(muHatLq, muStar);
    muHatLq = muHatLq*wLq;
    xHatLq = xHatLq*wLq;
    for i = 1:nBlock
        sigmaHatLq(:, :, i) = wLq'*squeeze(sigmaHatLq(:, :, i))*wLq;
    end
    
    % Rotate xHatMean0 to match muStar
    wMean0 = procrustes(muHatMean0, muStar);
    muHatMean0 = muHatMean0*wMean0;
    xHatMean0 = xHatMean0*wMean0;
    for i = 1:nBlock
        sigmaHatMean0(:, :, i) = wMean0'*squeeze(sigmaHatMean0(:, :, i))*wMean0;
    end
    
    % Rotate xHatLq0 to match muStar
    wLq0 = procrustes(muHatLq0, muStar);
    muHatLq0 = muHatLq0*wLq0;
    xHatLq0 = xHatLq0*wLq0;
    for i = 1:nBlock
        sigmaHatLq0(:, :, i) = wLq0'*squeeze(sigmaHatLq0(:, :, i))*wLq0;
    end
    
    
    
    
    
    scatterplot(xHatMean, tauHatMean, muHatMean, sigmaHatMean, 'r', 'Mean');
    hold on;
    scatterplot(xHatLq, tauHatLq, muHatLq, sigmaHatLq, 'b', 'Lq');
    
    scatterplot(xHatMean0, tauHatMean0, muHatMean0, sigmaHatMean0, 'g', 'Mean0');
    scatterplot(xHatLq0, tauHatLq0, muHatLq0, sigmaHatLq0, 'y', 'Lq0');
    
    plot(muStar(:, 1), muStar(:, 2), 'ko');
    
    plotLegend = legend('Mean: embedded points', ...
        'Mean: mean of cluster', 'Mean: 95% ellipse', ...
        'Lq: embedded points', 'Lq: mean of cluster', ...
        'Lq: 95% ellipse', 'Mean0: embedded points', ...
        'Mean0: mean of cluster', 'Mean0: 95% ellipse', ...
        'Lq0: embedded points', 'Lq0: mean of cluster', ...
        'Lq0: 95% ellipse', 'muStar');
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
