function ellipseplot(muHat, sigmaHat, color, hasLegend)

[eVec, eVal] = eig(sigmaHat);

[ind, ~] = find(eVal == max(max(eVal)));
eVecMax = eVec(:, ind);

phi = atan2(eVecMax(2), eVecMax(1));
if(phi < 0)
    phi = phi + 2*pi;
end

eValMax = max(max(eVal));
if (ind == 1)
    eValMin = eVal(2, 2);
else
    eValMin = eVal(1, 1);
end

a=2.4477*sqrt(eValMax);
b=2.4477*sqrt(eValMin);

gridTheta = linspace(0,2*pi);

xEllipse  = a*cos(gridTheta);
yEllipse  = b*sin(gridTheta);

R = [cos(phi) sin(phi); -sin(phi) cos(phi)];

rotateEllipse = [xEllipse;yEllipse]'*R;

hold on;
plot(rotateEllipse(:,1) + muHat(1), rotateEllipse(:,2) + muHat(2), ...
    strcat('-', color), 'HandleVisibility', hasLegend);

end