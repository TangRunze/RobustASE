function val = objectivefun(theta, q, xVec)
% Objective functions in optimization problem

val = sum((theta.^xVec./factorial(xVec)).^(1-q).*(xVec - theta));

end