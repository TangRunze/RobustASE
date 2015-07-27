
clear all
close all

%% --- Parameters Setting ---

seed = 12345;
rng(seed);

hasPlot = 1;

epsilon = 0.2;
muB = 5;
epsilonInB = 3;
scaleB = 1;
c = 0.1;
m = 20;
q = 0.8;

nVertex = 50;
nBlock = 2;
dimLatentPosition = nBlock;
rho = repmat(1/nBlock, 1, nBlock);
iStart = 1;
iEnd = 1;

if (hasPlot == 1)
    iEnd = 1;
end

% block probability matrix
B = (muB - epsilonInB)*ones(nBlock, nBlock) + 2*epsilonInB*eye(nBlock);
B = scaleB^2*B;
[U, S, V] = svd(B);
nuStar = U*sqrt(S);

% plot(muStar(:, 1), muStar(:, 2), '.b');
% hold on;
theta = -3*pi/4;
rotationMatrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];
nuStar = nuStar*rotationMatrix';
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
    PLq = zeros(nVertex, nVertex);
    PL1 = zeros(nVertex, nVertex);
    adjMatrixTotal0 = zeros(m, nVertex, nVertex);
    adjMatrixSum0 = zeros(nVertex, nVertex);
    PLq0 = zeros(nVertex, nVertex);
    PL10 = zeros(nVertex, nVertex);

    % Generate block assignment.
    tauStarInd = mnrnd(1, rho, nVertex);
    tauStar = (1:nBlock)*tauStarInd';
    
    for iGraph = ((iIter - 1)*m+1):(iIter*m)
        % Generate data if there does not exist one, otherwise read the
        % existing data.
        [adjMatrix, ~, adjMatrix0] = datagenerator_exp(nVertex, B, ...
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
            PLq(i, j) = fsolve(@(theta) ...
                objectivefun_exp(theta, q, squeeze(adjMatrixTotal(:,i,j))), ...
                1/adjMatrixMean(i, j));
            PLq(j, i) = PLq(i, j);
            PLq0(i, j) = fsolve(@(theta) ...
                objectivefun_exp(theta, q, squeeze(adjMatrixTotal0(:,i,j))), ...
                1/adjMatrixMean0(i, j));
            PLq0(j, i) = PLq0(i, j);
            PL1(i, j) = 1/adjMatrixMean(i, j);
            PL1(j, i) = PL1(i, j);
            PL10(i, j) = 1/adjMatrixMean0(i, j);
            PL10(j, i) = PL10(i, j);
        end
    end
    
    PStar = B(tauStar, tauStar);
    PStar(1:(nVertex+1):end) = 0;
    
    % Diagonal Augmentation.
    PL1DA = PL1;
    PL1DA = PL1DA + diag(sum(PL1DA, 2))/...
        (size(PL1DA, 2) - 1);
    
    PLqDA = PLq;
    PLqDA = PLqDA + diag(sum(PLqDA, 2))/...
        (size(PLqDA, 2) - 1);
    
    PL1DA0 = PL10;
    PL1DA0 = PL1DA0 + diag(sum(PL1DA0, 2))/...
        (size(PL1DA0, 2) - 1);
    
    PLqDA0 = PLq0;
    PLqDA0 = PLqDA0 + diag(sum(PLqDA0, 2))/...
        (size(PLqDA0, 2) - 1);

    
    
    % GMM o ASE
    xHatL1 = asge(PL1DA, dimLatentPosition);
    gmL1 = fitgmdist(xHatL1, nBlock, 'Replicates', 10);
    tauHatL1 = cluster(gmL1, xHatL1)';
    pTauHatL1 = posterior(gmL1, xHatL1)';
    muHatL1 = gmL1.mu;
    sigmaHatL1 = gmL1.Sigma;
    
    xHatLq = asge(PLqDA, dimLatentPosition);
    gmLq = fitgmdist(xHatLq, nBlock, 'Replicates', 10);
    tauHatLq = cluster(gmLq, xHatLq)';
    pTauHatLq = posterior(gmLq, xHatLq)';
    muHatLq = gmLq.mu;
    sigmaHatLq = gmLq.Sigma;
    
    
    xHatL10 = asge(PL1DA0, dimLatentPosition);
    gmL10 = fitgmdist(xHatL10, nBlock, 'Replicates', 10);
    tauHatL10 = cluster(gmL10, xHatL10)';
    pTauHatL10 = posterior(gmL10, xHatL10)';
    muHatL10 = gmL10.mu;
    sigmaHatL10 = gmL10.Sigma;
    
    xHatLq0 = asge(PLqDA0, dimLatentPosition);
    gmLq0 = fitgmdist(xHatLq0, nBlock, 'Replicates', 10);
    tauHatLq0 = cluster(gmLq0, xHatLq0)';
    pTauHatLq0 = posterior(gmLq0, xHatLq0)';
    muHatLq0 = gmLq0.mu;
    sigmaHatLq0 = gmLq0.Sigma;
    
    % Rotate xHatMean to match muStar
    wL1 = procrustes(muHatL1, nuStar);
    muHatL1 = muHatL1*wL1;
    xHatL1 = xHatL1*wL1;
    for i = 1:nBlock
        sigmaHatL1(:, :, i) = wL1'*squeeze(sigmaHatL1(:, :, i))*wL1;
    end
    
    % Rotate xHatLq to match muStar
    wLq = procrustes(muHatLq, nuStar);
    muHatLq = muHatLq*wLq;
    xHatLq = xHatLq*wLq;
    for i = 1:nBlock
        sigmaHatLq(:, :, i) = wLq'*squeeze(sigmaHatLq(:, :, i))*wLq;
    end
    
    % Rotate xHatMean0 to match muStar
    wL10 = procrustes(muHatL10, nuStar);
    muHatL10 = muHatL10*wL10;
    xHatL10 = xHatL10*wL10;
    for i = 1:nBlock
        sigmaHatL10(:, :, i) = wL10'*squeeze(sigmaHatL10(:, :, i))*wL10;
    end
    
    % Rotate xHatLq0 to match muStar
    wLq0 = procrustes(muHatLq0, nuStar);
    muHatLq0 = muHatLq0*wLq0;
    xHatLq0 = xHatLq0*wLq0;
    for i = 1:nBlock
        sigmaHatLq0(:, :, i) = wLq0'*squeeze(sigmaHatLq0(:, :, i))*wLq0;
    end
    
    
    % Calculate Rotation o GMM o ASE(expect value of MLqE)
    BL1Exp = [0.47619, 0.416667; 0.416667, 0.47619];
    PL1Exp = BL1Exp(tauStar, tauStar);
    PL1ExpDA = PL1Exp;
    PL1ExpDA = PL1ExpDA + diag(sum(PL1ExpDA, 2))/...
        (size(PL1ExpDA, 2) - 1);
    xHatL1Exp = asge(PL1Exp, dimLatentPosition);
    
    BLqExp = [9.27334, 1.9577; 1.9577, 9.27334];
    PLqExp = BLqExp(tauStar, tauStar);
    PLqExpDA = PLqExp;
    PLqExpDA = PLqExpDA + diag(sum(PLqExpDA, 2))/...
        (size(PLqExpDA, 2) - 1);
    xHatLqExp = asge(PLqExp, dimLatentPosition);
    
    indXHat = zeros(1, nBlock);
    for k = 1:nBlock
        indXHat(k) = find(tauStar == k, 1);
    end
    
    wL1Exp = procrustes(xHatL1Exp(indXHat, :), nuStar);
    xHatL1Exp = xHatL1Exp*wL1Exp;
    
    wLqExp = procrustes(xHatLqExp(indXHat, :), nuStar);
    xHatLqExp = xHatLqExp*wLqExp;
    
    
    % Plot
    if (hasPlot == 1)
        scatterplot(xHatL1, tauHatL1, muHatL1, sigmaHatL1, 'r', 'Mean');
        hold on;
        scatterplot(xHatLq, tauHatLq, muHatLq, sigmaHatLq, 'b', 'Lq');
        scatterplot(xHatL10, tauHatL10, muHatL10, sigmaHatL10, 'g', 'Mean0');
        scatterplot(xHatLq0, tauHatLq0, muHatLq0, sigmaHatLq0, 'm', 'Lq0');
        plot(xHatL1Exp(:, 1), xHatL1Exp(:, 2), 'oc', 'MarkerEdgeColor', 'c',...
            'MarkerFaceColor', 'c', 'markersize', 6);
        plot(xHatLqExp(:, 1), xHatLqExp(:, 2), 'oy', 'MarkerEdgeColor', 'y',...
            'MarkerFaceColor', 'y', 'markersize', 6);
        
        plot(nuStar(:, 1), nuStar(:, 2), 'ks', 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'k', 'markersize', 5);
        
        plotLegend = legend(['$\hat{X}_{q=1}$ \ \ with $\epsilon=' num2str(epsilon) '$'], ...
            ['$\hat{\nu}_{q=1}$ \ \ with $\epsilon=' num2str(epsilon) '$'], ...
            ['$\hat{\Sigma}_{q=1}$ \ \ with $\epsilon=' num2str(epsilon) '$'], ...
            ['$\hat{X}_{q=' num2str(q) '}$ \ \ with $\epsilon=' num2str(epsilon) '$'], ...
            ['$\hat{\nu}_{q=' num2str(q) '}$ \ \ with $\epsilon=' num2str(epsilon) '$'], ...
            ['$\hat{\Sigma}_{q=' num2str(q) '}$ \ \ with $\epsilon=' num2str(epsilon) '$'], ...
            '$\hat{X}_{q=1}$ \ \ with $\epsilon=0$', ...
            '$\hat{\nu}_{q=1}$ \ \ with $\epsilon=0$', ...
            '$\hat{\Sigma}_{q=1}$ \ \ with $\epsilon=0$', ...
            ['$\hat{X}_{q=' num2str(q) '}$ \ \ with $\epsilon=0$'], ...
            ['$\hat{\nu}_{q=' num2str(q) '}$ \ \ with $\epsilon=0$'], ...
            ['$\hat{\Sigma}_{q=' num2str(q) '}$ \ \ with $\epsilon=0$'], ...
            ['$\nu_{q=1}^{\mathrm{Expect}}$ \ \ with $\epsilon=' num2str(epsilon) '$'], ...
            ['$\nu_{q=' num2str(q) '}^{\mathrm{Expect}}$ \ \ with $\epsilon=' num2str(epsilon) '$'], ...
            '$\nu$');
        set(plotLegend, 'FontSize', 14, 'Interpreter', 'latex');
        hold off;
        plotTitle = title(strcat('n=', num2str(nVertex), ', K=', ...
            num2str(nBlock), ', $\epsilon$=', num2str(epsilon), ', diag=', ...
            num2str(muB + epsilonInB), ', offdiag=', ...
            num2str(muB - epsilonInB), ', q=', num2str(q), ...
            ', C=', num2str(c), ', m=', num2str(m), ...
            ', seed=', num2str(seed)));
        set(plotTitle, 'FontSize', 14, 'Interpreter', 'latex');
    end
    
    % mean
    ADiffMean(iIter) = norm(adjMatrixMean - PStar);
    
%     xHatMean = asge(adjMatrixMean, dimLatentPosition);
%     gm = fitgmdist(xHatMean, nBlock, 'Replicates', 10);
%     tauHatMean = cluster(gm, xHatMean)';
    errorRateMean(iIter) = errorratecalculator(tauStar, tauHatL1, nBlock);
    
    AXDiffMean(iIter) = norm(xHatL1*xHatL1' - PStar);
    
    % Lq
    ADiffLq(iIter) = norm(PLq - PStar);
    
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
