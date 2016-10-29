function val = mlqe_obj_lognormal_root(theta, xVec, q)

mu = theta(1)^2;
sigma = theta(2)^2;

fVal = 1/sigma/sqrt(2*pi)./xVec.*exp(-(log(xVec) - mu).^2/2/(sigma^2));
fqVal = fVal.^(1-q);

val = sum(fqVal.*(log(xVec) - mu)) + ...
    sum(fqVal.*(sigma^2 - (log(xVec) - mu).^2))^2;