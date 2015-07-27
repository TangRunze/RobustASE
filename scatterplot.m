function scatterplot(xHat, tauHat, muHat, sigmaHat, color, comment)

    plot(xHat(:, 1), xHat(:, 2), strcat('.', color));
    hold on;
    plot(muHat(:, 1), muHat(:, 2), strcat('s', color), 'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', color, 'MarkerSize', 5);
    
    ellipseplot(squeeze(muHat(1, :)), squeeze(sigmaHat(:, :, 1)), color,...
        'on');

    ellipseplot(squeeze(muHat(2, :)), squeeze(sigmaHat(:, :, 2)), color,...
        'off');
end