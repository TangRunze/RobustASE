function [adjMatrix, tauStar, vertexNoiseInd] = ...
    datagenerator_conti(nVertex, nBlock, B, rho, epsilon, scaleVar, iGraph)

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
    vertexNoise = (binornd(1, epsilon, 1, nVertex) == 1);
    vertexNoiseInd = find(vertexNoise);
    
    muMatrix = B(tauStar, tauStar);
    sigmaMatrix = muMatrix/scaleVar;
    sigmaMatrix(vertexNoise, :) = ...
        sigmaMatrix(vertexNoise, :)*scaleVar;
    sigmaMatrix(:, vertexNoise) = ...
        sigmaMatrix(:, vertexNoise)*scaleVar;
    sigmaMatrix(vertexNoise, vertexNoise) = ...
        sigmaMatrix(vertexNoise, vertexNoise)/scaleVar;

    % Generate graph adjacency matrix.
    adjMatrix = zeros(nVertex, nVertex);
    for i = 1:nVertex
        for j = (i+1):nVertex
            adjMatrix(i, j) = normrnd(muMatrix(i, j), ...
                sigmaMatrix(i, j));
            while (adjMatrix(i, j) < 0) || ...
                    (adjMatrix(i, j) > 2*muMatrix(i, j))
                adjMatrix(i, j) = normrnd(muMatrix(i, j), ...
                    sigmaMatrix(i, j));
            end
            adjMatrix(j, i) = adjMatrix(i, j);
        end
    end
    
    % Save the data
    save(['data/sim-n' num2str(nVertex) '-diag' num2str(B(1, 1)) ...
        '-offdiag' num2str(B(1, 2)) '-eps' num2str(epsilon) ...
        '-graph' int2str(iGraph) '.mat'], 'adjMatrix', 'tauStar', ...
        'vertexNoiseInd');
else
    % Read the existing data
    data = load(['data/sim-n' num2str(nVertex) '-diag' num2str(B(1, 1)) ...
        '-offdiag' num2str(B(1, 2)) '-eps' num2str(epsilon) ...
        '-graph' int2str(iGraph) '.mat']);
    adjMatrix = data.adjMatrix;
    tauStar = data.tauStar;
    vertexNoiseInd = data.vertexNoiseInd;
end
