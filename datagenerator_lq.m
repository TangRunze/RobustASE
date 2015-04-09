function [adjMatrix, edgeNoise, adjMatrix0] = ...
    datagenerator_lq(nVertex, B, tauStar, epsilon, c, iGraph)

% Generate data if there does not exist one, otherwise read the
% existing data.
if exist(['data/sim-n' num2str(nVertex) '-diag' num2str(B(1, 1)) ...
        '-offdiag' num2str(B(1, 2)) '-eps' num2str(epsilon) ...
        '-graph' int2str(iGraph) '.mat'], 'file') == 0
    
    disp(['Generating graph ' int2str(iGraph) '...'])
    
    % Edge with noise.
    edgeNoise = (binornd(1, epsilon, nVertex, nVertex) == 1);
    edgeNoise = triu(edgeNoise, 1);
    edgeNoise = edgeNoise + edgeNoise';
    
    % Generate graph adjacency matrix.
    adjMatrix = zeros(nVertex, nVertex);
    adjMatrix0 = adjMatrix;
    for i = 1:nVertex
        for j = (i+1):nVertex
            mu = B(tauStar(i), tauStar(j));
            adjMatrix0(i, j) = poissrnd(mu);
            if (edgeNoise(i, j) == 1)
                mu = c*rand;
                adjMatrix(i, j) = poissrnd(mu);
            else
                adjMatrix(i, j) = adjMatrix0(i, j);
            end
            adjMatrix(j, i) = adjMatrix(i, j);
            adjMatrix0(j, i) = adjMatrix0(i, j);
        end
    end
    
    % Save the data
    save(['data/sim-n' num2str(nVertex) '-diag' num2str(B(1, 1)) ...
        '-offdiag' num2str(B(1, 2)) '-eps' num2str(epsilon) ...
        '-graph' int2str(iGraph) '.mat'], 'adjMatrix', 'tauStar', ...
        'edgeNoise', 'adjMatrix0');
else
    % Read the existing data
    data = load(['data/sim-n' num2str(nVertex) '-diag' num2str(B(1, 1)) ...
        '-offdiag' num2str(B(1, 2)) '-eps' num2str(epsilon) ...
        '-graph' int2str(iGraph) '.mat']);
    adjMatrix = data.adjMatrix;
    % tauStar = data.tauStar;
    edgeNoise = data.edgeNoise;
    adjMatrix0 = data.adjMatrix0;
end
