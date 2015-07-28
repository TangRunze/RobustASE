
clear all
close all

%% --- Parameters Setting ---

seed = 12345;
rng(seed);

hasPlot = 1;
nGraph = 100;

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


if (hasPlot == 1)
    nGraph = 1;
end

% block probability matrix
B = (muB - epsilonInB)*ones(nBlock, nBlock) + 2*epsilonInB*eye(nBlock);
B = scaleB^2*B;
[U, S, V] = svd(B);
nuStar = U*sqrt(S);

theta = -3*pi/4;
rotationMatrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];
nuStar = nuStar*rotationMatrix';

%% --- Parallel Computing ---
% delete(gcp('nocreate'))
% parpool(nCore);


%% Monte Carlo Simulation

errorRateL1 = zeros(1, nGraph);
errorRateLq = zeros(1, nGraph);
PDiffL1 = zeros(1, nGraph);
PDiffLq = zeros(1, nGraph);
XXTDiffL1 = zeros(1, nGraph);
XXTDiffLq = zeros(1, nGraph);
nunuTDiffL1 = zeros(1, nGraph);
nunuTDiffLq = zeros(1, nGraph);
errorRateL10 = zeros(1, nGraph);
errorRateLq0 = zeros(1, nGraph);
PDiffL10 = zeros(1, nGraph);
PDiffLq0 = zeros(1, nGraph);
XXTDiffL10 = zeros(1, nGraph);
XXTDiffLq0 = zeros(1, nGraph);
nunuTDiffL10 = zeros(1, nGraph);
nunuTDiffLq0 = zeros(1, nGraph);

for iIter = 1:nGraph
    iIter
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
        [adjMatrix, ~, adjMatrix0, tauStar] = datagenerator_exp(nVertex,...
            B, tauStar, epsilon, c, iGraph);
        adjMatrixTotal(iGraph - (iIter - 1)*m, :, :) = adjMatrix;
        adjMatrixSum = adjMatrixSum + adjMatrix;
        adjMatrixTotal0(iGraph - (iIter - 1)*m, :, :) = adjMatrix0;
        adjMatrixSum0 = adjMatrixSum0 + adjMatrix0;
    end
    
    adjMatrixMean = adjMatrixSum/m;
    adjMatrixMean0 = adjMatrixSum0/m;
    for i = 1:nVertex
        for j = (i + 1):nVertex
%             [i j]
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
    nuHatL1 = gmL1.mu;
    sigmaHatL1 = gmL1.Sigma;
    
    xHatLq = asge(PLqDA, dimLatentPosition);
    gmLq = fitgmdist(xHatLq, nBlock, 'Replicates', 10);
    tauHatLq = cluster(gmLq, xHatLq)';
    pTauHatLq = posterior(gmLq, xHatLq)';
    nuHatLq = gmLq.mu;
    sigmaHatLq = gmLq.Sigma;
    
    
    xHatL10 = asge(PL1DA0, dimLatentPosition);
    gmL10 = fitgmdist(xHatL10, nBlock, 'Replicates', 10);
    tauHatL10 = cluster(gmL10, xHatL10)';
    pTauHatL10 = posterior(gmL10, xHatL10)';
    nuHatL10 = gmL10.mu;
    sigmaHatL10 = gmL10.Sigma;
    
    xHatLq0 = asge(PLqDA0, dimLatentPosition);
    gmLq0 = fitgmdist(xHatLq0, nBlock, 'Replicates', 10);
    tauHatLq0 = cluster(gmLq0, xHatLq0)';
    pTauHatLq0 = posterior(gmLq0, xHatLq0)';
    nuHatLq0 = gmLq0.mu;
    sigmaHatLq0 = gmLq0.Sigma;
    
    % Plot
    if (hasPlot == 1)
        
        % Rotate xHatL1 to match muStar
        wL1 = procrustes(nuHatL1, nuStar);
        nuHatL1 = nuHatL1*wL1;
        xHatL1 = xHatL1*wL1;
        for i = 1:nBlock
            sigmaHatL1(:, :, i) = wL1'*squeeze(sigmaHatL1(:, :, i))*wL1;
        end
        
        % Rotate xHatLq to match muStar
        wLq = procrustes(nuHatLq, nuStar);
        nuHatLq = nuHatLq*wLq;
        xHatLq = xHatLq*wLq;
        for i = 1:nBlock
            sigmaHatLq(:, :, i) = wLq'*squeeze(sigmaHatLq(:, :, i))*wLq;
        end
        
        % Rotate xHatL10 to match muStar
        wL10 = procrustes(nuHatL10, nuStar);
        nuHatL10 = nuHatL10*wL10;
        xHatL10 = xHatL10*wL10;
        for i = 1:nBlock
            sigmaHatL10(:, :, i) = wL10'*squeeze(sigmaHatL10(:, :, i))*wL10;
        end
        
        % Rotate xHatLq0 to match muStar
        wLq0 = procrustes(nuHatLq0, nuStar);
        nuHatLq0 = nuHatLq0*wLq0;
        xHatLq0 = xHatLq0*wLq0;
        for i = 1:nBlock
            sigmaHatLq0(:, :, i) = wLq0'*squeeze(sigmaHatLq0(:, :, i))*wLq0;
        end
        
        
        % Calculate Rotation o GMM o ASE(expect value of MLqE)
        BL1Exp = arrayfun(@(x) l1expectation_exp(x, c, epsilon), B);
        PL1Exp = BL1Exp(tauStar, tauStar);
        PL1ExpDA = PL1Exp;
        PL1ExpDA(1:(nVertex+1):end) = 0;
        PL1ExpDA = PL1ExpDA + diag(sum(PL1ExpDA, 2))/...
            (size(PL1ExpDA, 2) - 1);
        xHatL1Exp = asge(PL1Exp, dimLatentPosition);
        
        BLqExp = arrayfun(@(x) lqexpectation_exp(x, q, c, epsilon), B);
        PLqExp = BLqExp(tauStar, tauStar);
        PLqExpDA = PLqExp;
        PLqExpDA(1:(nVertex+1):end) = 0;
        PLqExpDA = PLqExpDA + diag(sum(PLqExpDA, 2))/...
            (size(PLqExpDA, 2) - 1);
        xHatLqExp = asge(PLqExp, dimLatentPosition);
        
        BL1Exp0 = arrayfun(@(x) l1expectation_exp(x, c, 0), B);
        PL1Exp0 = BL1Exp0(tauStar, tauStar);
        PL1ExpDA0 = PL1Exp0;
        PL1ExpDA0(1:(nVertex+1):end) = 0;
        PL1ExpDA0 = PL1ExpDA0 + diag(sum(PL1ExpDA0, 2))/...
            (size(PL1ExpDA0, 2) - 1);
        xHatL1Exp0 = asge(PL1Exp0, dimLatentPosition);
        
        BLqExp0 = arrayfun(@(x) lqexpectation_exp(x, q, c, 0), B);
        PLqExp0 = BLqExp0(tauStar, tauStar);
        PLqExpDA0 = PLqExp0;
        PLqExpDA0(1:(nVertex+1):end) = 0;
        PLqExpDA0 = PLqExpDA0 + diag(sum(PLqExpDA0, 2))/...
            (size(PLqExpDA0, 2) - 1);
        xHatLqExp0 = asge(PLqExp0, dimLatentPosition);
        
        
        indXHat = zeros(1, nBlock);
        for k = 1:nBlock
            indXHat(k) = find(tauStar == k, 1);
        end
        
        wL1Exp = procrustes(xHatL1Exp(indXHat, :), nuStar);
        xHatL1Exp = xHatL1Exp*wL1Exp;
        
        wLqExp = procrustes(xHatLqExp(indXHat, :), nuStar);
        xHatLqExp = xHatLqExp*wLqExp;
        
        wL1Exp0 = procrustes(xHatL1Exp0(indXHat, :), nuStar);
        xHatL1Exp0 = xHatL1Exp0*wL1Exp0;
        
        wLqExp0 = procrustes(xHatLqExp0(indXHat, :), nuStar);
        xHatLqExp0 = xHatLqExp0*wLqExp0;
        
        % Plot
        plot(nuStar(:, 1), nuStar(:, 2), 'ks', 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'k', 'markersize', 5);
        hold on;
        scatterplot(xHatL1, tauHatL1, nuHatL1, sigmaHatL1, 'r', 'Mean');

        scatterplot(xHatLq, tauHatLq, nuHatLq, sigmaHatLq, 'b', 'Lq');
        scatterplot(xHatL10, tauHatL10, nuHatL10, sigmaHatL10, 'g', 'Mean0');
        scatterplot(xHatLq0, tauHatLq0, nuHatLq0, sigmaHatLq0, 'm', 'Lq0');
        plot(xHatL1Exp(:, 1), xHatL1Exp(:, 2), '^c', ...
            'markersize', 6);
        plot(xHatLqExp(:, 1), xHatLqExp(:, 2), '^k', ...
            'markersize', 6);
        plot(xHatL1Exp0(:, 1), xHatL1Exp0(:, 2), 'xc', ...
            'markersize', 6);
        plot(xHatLqExp0(:, 1), xHatLqExp0(:, 2), 'xk', ...
            'markersize', 6);
        

        
        plotLegend = legend('$\nu$', ...
            ['$\hat{X}_{q=1}$ \ \ with $\epsilon=' num2str(epsilon) '$'], ...
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
            ['$\nu_{q=1}^{\mathrm{Expect}}$ \ \ with $\epsilon=0$'], ...
            ['$\nu_{q=' num2str(q) '}^{\mathrm{Expect}}$ \ \ with $\epsilon=0$']);
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
    
    % L1
    PDiffL1(iIter) = norm(PL1 - PStar, 'fro');
    errorRateL1(iIter) = errorratecalculator(tauStar, tauHatL1, nBlock);
    XXTDiffL1(iIter) = norm(xHatL1*xHatL1' - PStar, 'fro');
    nunuTDiffL1(iIter) = norm(nuHatL1(tauHatL1, :)*nuHatL1(tauHatL1, :)'...
        - PStar, 'fro');
    
    % Lq
    PDiffLq(iIter) = norm(PLq - PStar, 'fro');
    errorRateLq(iIter) = errorratecalculator(tauStar, tauHatLq, nBlock);
    XXTDiffLq(iIter) = norm(xHatLq*xHatLq' - PStar, 'fro');
    nunuTDiffLq(iIter) = norm(nuHatLq(tauHatLq, :)*nuHatLq(tauHatLq, :)'...
        - PStar, 'fro');
    
    % L10
    PDiffL10(iIter) = norm(PL10 - PStar, 'fro');
    errorRateL10(iIter) = errorratecalculator(tauStar, tauHatL10, nBlock);
    XXTDiffL10(iIter) = norm(xHatL10*xHatL10' - PStar, 'fro');
    nunuTDiffL10(iIter) = norm(nuHatL10(tauHatL10, :)*...
        nuHatL10(tauHatL10, :)' - PStar, 'fro');
    
    % Lq0
    PDiffLq0(iIter) = norm(PLq0 - PStar, 'fro');
    errorRateLq0(iIter) = errorratecalculator(tauStar, tauHatLq0, nBlock);
    XXTDiffLq0(iIter) = norm(xHatLq0*xHatLq0' - PStar, 'fro');
    nunuTDiffLq0(iIter) = norm(nuHatLq0(tauHatLq0, :)*...
        nuHatLq0(tauHatLq0, :)' - PStar, 'fro');
end

[mean(errorRateL1), mean(errorRateLq)]
pValueError = signtest(errorRateL1, errorRateLq)

[mean(PDiffL1), mean(PDiffLq)]
pValuePDiff = signtest(PDiffL1, PDiffLq)

[mean(XXTDiffL1), mean(XXTDiffLq)]
pValueXXTDiff = signtest(XXTDiffL1, XXTDiffLq)

[mean(nunuTDiffL1), mean(nunuTDiffLq)]
pValuenunuTDiff = signtest(nunuTDiffL1, nunuTDiffLq)

signtest(PDiffLq, XXTDiffLq)
signtest(XXTDiffLq, nunuTDiffLq)



[mean(errorRateL10), mean(errorRateLq0)]
pValueError0 = signtest(errorRateLq0, errorRateL10)

[mean(PDiffL10), mean(PDiffLq0)]
pValuePDiff0 = signtest(PDiffLq0, PDiffL10)

[mean(XXTDiffL10), mean(XXTDiffLq0)]
pValueXXTDiff0 = signtest(XXTDiffLq0, XXTDiffL10)

[mean(nunuTDiffL10), mean(nunuTDiffLq0)]
pValuenunuTDiff0 = signtest(nunuTDiffLq0, nunuTDiffL10)

signtest(PDiffLq0, XXTDiffLq0)
signtest(XXTDiffLq0, nunuTDiffLq0)

