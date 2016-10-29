

dataName = 'CPAC200';
qVec = [1, 0.95, 0.9, 0.8, 0.7];
% qVec = [1, 0.95, 0.9, 0.85];

pairVec = [3, 4];
inputFileName = ['../../Data/', dataName, '_pair_', num2str(pairVec(1)),...
    '_', num2str(pairVec(2)), '.csv'];

data = csvread(inputFileName);

% nv = ((data == 0) | (data == 1));
% data(nv) = [];
data = data + 2;
data = log(data);

histogram(data, 'Normalization', 'pdf');

muVec = [];
sigmaVec = [];
xVec = (1:(2*max(data)*100))/100;

% options = optimoptions('fminunc', 'TolX', 1e-6, ...
%     'MaxIter', 10000, 'MaxFunEvals', 10000);

muHat = mean(log(data));
sigmaHat = sqrt(mean((log(data) - muHat).^2));

thetaInit = [muHat, sqrt(sigmaHat)];

for q = qVec
    [thetaHat, fVal] = fminunc(@(theta) ...
        mlqe_obj_lognormal(theta, data, q), thetaInit);
    muHat = thetaHat(1);
    sigmaHat = thetaHat(2)^2;

    muVec = [muVec, muHat];
    sigmaVec = [sigmaVec, sigmaHat];
    
    % histfit(data)
    hold on;
    fHat = 1/sigmaHat/sqrt(2*pi)./xVec.*exp(-(log(xVec) - muHat).^2/2/(sigmaHat^2));
    plot(xVec, fHat);
end

% num2str(qVec', muVec', sigmaVec, 'q=%.2f, mu=%.2f, sigma=%.2f')
% [num2str(qVec', 'q=%.2f'), num2str(muVec', 'mu=%.2f')]
legendCell = cellstr(num2str(qVec', 'q=%.2f'));
legendCell = {'Raw Data', legendCell{1:length(legendCell)}};
for i = 1:(length(legendCell) - 1)
    legendCell{i+1} = [legendCell{i+1}, ', mu=', ...
        num2str(muVec(i), '%.2f'), ', sigma=', num2str(sigmaVec(i), '%.2f')];
end
legend(legendCell, 'Location', 'SouthOutside');
title(['MLqE for Element (', num2str(pairVec(1)), ', ', num2str(pairVec(2)), ')']);

figureFileName = [dataName, '_pair_', num2str(pairVec(1)),...
    '_', num2str(pairVec(2)), '.png'];

set(gcf, 'PaperUnits', 'inches', 'PaperPosition',[0 0 15 9])
print('-dpng', figureFileName, '-r100');
