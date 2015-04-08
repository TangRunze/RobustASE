function [adjMatrix, edgeNoise] = ...
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
    for i = 1:nVertex
        for j = (i+1):nVertex
            if (edgeNoise(i, j) == 1)
                mu = c*rand;
            else
                mu = B(tauStar(i), tauStar(j));
            end
            adjMatrix(i, j) = poissrnd(mu);
            adjMatrix(j, i) = adjMatrix(i, j);
        end
    end
    
    % Save the data
    save(['data/sim-n' num2str(nVertex) '-diag' num2str(B(1, 1)) ...
        '-offdiag' num2str(B(1, 2)) '-eps' num2str(epsilon) ...
        '-graph' int2str(iGraph) '.mat'], 'adjMatrix', 'tauStar', ...
        'edgeNoise');
else
    % Read the existing data
    data = load(['data/sim-n' num2str(nVertex) '-diag' num2str(B(1, 1)) ...
        '-offdiag' num2str(B(1, 2)) '-eps' num2str(epsilon) ...
        '-graph' int2str(iGraph) '.mat']);
    adjMatrix = data.adjMatrix;
    % tauStar = data.tauStar;
    edgeNoise = data.edgeNoise;
end
