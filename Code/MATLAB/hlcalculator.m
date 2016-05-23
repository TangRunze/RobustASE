function adjMatrixHL = hlcalculator(adjMatrixTotal, m)

adjMatrixPair = zeros(m*(m+1)/2, size(adjMatrixTotal, 2), ...
    size(adjMatrixTotal, 3));

k = 0;
for i = 1:m
    for j = i:m
        k = k + 1;
        adjMatrixPair(k, :, :) = (squeeze(adjMatrixTotal(i, :, :)) + ...
            squeeze(adjMatrixTotal(j, :, :)))/2;
    end
end

adjMatrixHL = squeeze(median(adjMatrixPair));

% Version 2
% adjMatrix = zeros(m, size(adjMatrixTmp, 1), size(adjMatrixTmp, 1));
% adjMatrixAvg = zeros(m*(m+1)/2, size(adjMatrixTmp, 1), size(adjMatrixTmp, 1));
% 
% for i = 1:m
%     adjMatrix(i, :, :) = (adjMatrixTmp >= i);
% end
% 
% k = 0;
% for i = 1:m
%     for j = i:m
%         k = k + 1;
%         adjMatrixAvg(k, :, :) = (adjMatrix(i, :, :) + adjMatrix(j, :, :))/2;
%     end
% end
% 
% adjMatrixHL = squeeze(median(adjMatrixAvg));

% Version 1
% if mod(m*(m+1)/2, 2) == 0
%     nvl1 = adjMatrixTmp.*(adjMatrixTmp+1)/2 >= m*(m+1)/4 + 1;
%     nvr1 = adjMatrixTmp.*(adjMatrixTmp+1)/2 >= m*(m+1)/4;
%     nvl0 = (m-adjMatrixTmp).*(m-adjMatrixTmp+1)/2 >= m*(m+1)/4;
%     nvr0 = (m-adjMatrixTmp).*(m-adjMatrixTmp+1)/2 >= m*(m+1)/4 + 1;
%     adjMatrixHLleft = zeros(size(adjMatrixTmp));
%     adjMatrixHLleft(nvl1) = 1;
%     adjMatrixHLleft(nvl0) = 0;
%     adjMatrixHLleft((~nvl1) & (~nvl0)) = 0.5;
%     adjMatrixHLright = zeros(size(adjMatrixTmp));
%     adjMatrixHLright(nvr1) = 1;
%     adjMatrixHLright(nvr0) = 0;
%     adjMatrixHLright((~nvr1) & (~nvr0)) = 0.5;
%     adjMatrixHL = (adjMatrixHLleft + adjMatrixHLright)/2;
% else
%     nv1 = adjMatrixTmp.*(adjMatrixTmp+1)/2 >= (m*(m+1)+2)/4;
%     nv0 = (m-adjMatrixTmp).*(m-adjMatrixTmp+1)/2 >= (m*(m+1)+2)/4;
%     adjMatrixHL = zeros(size(adjMatrixTmp));
%     adjMatrixHL(nv1) = 1;
%     adjMatrixHL((~nv0) & (~nv1)) = 0.5;
% end