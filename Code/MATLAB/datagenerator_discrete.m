function [adjMatrix, tauStar, xStar, vertexNoise] = ...
    datagenerator_discrete(nVertex, nBlock, dimLatentPosition, B, nuStar, rho, ...
    epsilon, scale, ratio, iGraph)

% Generate data if there does not exist one, otherwise read the
% existing data.
if exist(['data/sim-n' num2str(nVertex) '-diag' num2str(B(1, 1)) ...
        '-offdiag' num2str(B(1, 2)) '-eps' num2str(epsilon) ...
        '-graph' int2str(iGraph) '.mat'], 'file') == 0
    
    disp(['Generating graph ' int2str(iGraph) '...'])
    
    % Assign block memberships.
    nVectorStar = rho*nVertex;
    
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
    
    % Vertices with noise.
    vertexNoise = find(binornd(1, epsilon, 1, nVertex) == 1);
    
    xStar = zeros(nVertex, dimLatentPosition);
    xStar(1:nVectorStarEnd(nBlock), :) = ...
        nuStar(tauStar(1:nVectorStarEnd(nBlock)), :);
    
%     % Generate true latent positions.
%     xNoise = xStar;
%     for iVertex = vertexNoise
%         xTmp = 2*rand(1, dimLatentPosition) - 1;
%         pXTmp = xTmp*xTmp';
%         while (pXTmp > 1) || (pXTmp < 0)
%             xTmp = 2*rand(1, dimLatentPosition) - 1;
%             pXTmp = xTmp*xTmp';
%         end
%         pTmp = xNoise*xTmp';
%         while any(pTmp > 1) || any(pTmp < 0)
%             xTmp = 2*rand(1, dimLatentPosition) - 1;
%             pXTmp = xTmp*xTmp';
%             while (pXTmp > 1) || (pXTmp < 0)
%                 xTmp = 2*rand(1, dimLatentPosition) - 1;
%                 pXTmp = xTmp*xTmp';
%             end
%             pTmp = xNoise*xTmp';
%         end
%         xNoise(iVertex, :) = xTmp;
%     end
    
    xNoise = xStar;
    for iVertex = vertexNoise
        xTmp = (2*rand(1, dimLatentPosition) - 1)*scale*ratio;
        pXTmp = xTmp*xTmp';
        while (pXTmp < 0)
            xTmp = (2*rand(1, dimLatentPosition) - 1)*scale*ratio;
            pXTmp = xTmp*xTmp';
        end
        pTmp = xNoise*xTmp';
        while any(pTmp < 0)
            xTmp = (2*rand(1, dimLatentPosition) - 1)*scale*ratio;
            pXTmp = xTmp*xTmp';
            while (pXTmp < 0)
                xTmp = (2*rand(1, dimLatentPosition) - 1)*scale*ratio;
                pXTmp = xTmp*xTmp';
            end
            pTmp = xNoise*xTmp';
        end
        xNoise(iVertex, :) = xTmp;
    end

    % Generate graph adjacency matrix.
    pMatrix = xNoise*xNoise';
%     adjMatrix = reshape(binornd(ones(1, nVertex*nVertex), ...
%         reshape(pMatrix, 1, nVertex*nVertex)), nVertex, nVertex);
    adjMatrix = poissrnd(pMatrix);
    adjMatrix = triu(adjMatrix, 1);
    adjMatrix = adjMatrix + adjMatrix';
    
    % Save the data
    save(['data/sim-n' num2str(nVertex) '-diag' num2str(B(1, 1)) ...
        '-offdiag' num2str(B(1, 2)) '-eps' num2str(epsilon) ...
        '-graph' int2str(iGraph) '.mat'], 'adjMatrix', 'tauStar', ...
        'xStar', 'vertexNoise');
else
    % Read the existing data
    data = load(['data/sim-n' num2str(nVertex) '-diag' num2str(B(1, 1)) ...
        '-offdiag' num2str(B(1, 2)) '-eps' num2str(epsilon) ...
        '-graph' int2str(iGraph) '.mat']);
    adjMatrix = data.adjMatrix;
    tauStar = data.tauStar;
    xStar = data.xStar;
    vertexNoise = data.vertexNoise;
end
