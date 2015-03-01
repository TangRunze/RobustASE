function [adjMatrixSum, tauStar, xStar] = datagenerator(nVertex, nBlock, ...
    dimLatentPosition, B, nuStar, rho, epsilon, m, iIter)
% Generate data if there does not exist one, otherwise read the
% existing data.

if exist(['data/sim-n' num2str(nVertex) '-diag' num2str(B(1, 1)) ...
        '-offdiag' num2str(B(1, 2)) '-eps' num2str(epsilon) ...
        '-size' num2str(m) '-group' int2str(iIter) '.mat'], 'file') == 0
    
    disp(['Generating graph group ' int2str(iIter) '...'])
    
    adjMatrixSum = zeros(nVertex, nVertex);
    
    % Assign block memberships randomly.
    nVectorStar = mnrnd(nVertex, [(1-epsilon)*rho, epsilon]);
    
    % Calculate the sizes
    nVectorStarStart = cumsum(nVectorStar);
    nVectorStarStart = [1, nVectorStarStart(1:(end-1)) + 1];
    nVectorStarEnd = cumsum(nVectorStar);
    
    % True block assignment tauStar (1-by-nVertex)
    tauStar = zeros(1, nVertex);
    for i = 1:nBlock
        tauStar(nVectorStarStart(i):nVectorStarEnd(i)) = ...
            i*ones(1, nVectorStar(i));
    end
    
    for iGraph = 1:m
        % Generate true latent positions.
        xStar = zeros(nVertex, dimLatentPosition);
        xStar(1:nVectorStarEnd(nBlock), :) = ...
            nuStar(tauStar(1:nVectorStarEnd(nBlock)), :);
        for iVertex = nVectorStarStart(nBlock+1):nVertex
            xTmp = 2*rand(1, dimLatentPosition) - 1;
            pXTmp = xTmp*xTmp';
            while (pXTmp > 1) || (pXTmp < 0)
                xTmp = rand(1, dimLatentPosition);
                pXTmp = xTmp*xTmp';
            end
            pTmp = xStar(1:(iVertex-1), :)*xTmp';
            while any(pTmp > 1) || any(pTmp < 0)
                xTmp = 2*rand(1, dimLatentPosition) - 1;
                pXTmp = xTmp*xTmp';
                while (pXTmp > 1) || (pXTmp < 0)
                    xTmp = rand(1, dimLatentPosition);
                    pXTmp = xTmp*xTmp';
                end
                pTmp = xStar(1:(iVertex-1), :)*xTmp';
            end
            xStar(iVertex, :) = xTmp;
        end
        
        % Generate graph adjacency matrix.
        pMatrix = xStar*xStar';
        adjMatrixTmp = reshape(binornd(ones(1, nVertex*nVertex), ...
            reshape(pMatrix, 1, nVertex*nVertex)), nVertex, nVertex);
        adjMatrixTmp = triu(adjMatrixTmp, 1);
        adjMatrixTmp = adjMatrixTmp + adjMatrixTmp';
        
        adjMatrixSum = adjMatrixSum + adjMatrixTmp;
        
    end
    
    % Save the data
    save(['data/sim-n' num2str(nVertex) '-diag' num2str(B(1, 1)) ...
        '-offdiag' num2str(B(1, 2)) '-eps' num2str(epsilon) ...
        '-size' num2str(m) '-group' int2str(iIter) '.mat'], ...
        'adjMatrixSum', 'tauStar');
else
    % Read the existing data
    data = load(['data/sim-n' num2str(nVertex) '-diag' num2str(B(1, 1)) ...
        '-offdiag' num2str(B(1, 2)) '-eps' num2str(epsilon) ...
        '-size' num2str(m) '-group' int2str(iIter) '.mat']);
    adjMatrixSum = data.adjMatrixSum;
    tauStar = data.tauStar;
end
