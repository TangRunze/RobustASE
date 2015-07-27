function val = objectivefun_exp(theta, q, xVec)
% Objective functions in optimization problem

val = sum((exp(-theta*xVec)*theta).^(1 - q).*(1/theta - xVec));
end