function theta = ltsolve(xVec, q)
    theta = mean(xVec);
    thetaMin = min(xVec);
    thetaMax = max(xVec);
    xFactorial = factorial(xVec);
    f = sum((q.^(exp(-theta)*theta.^xVec./xFactorial)).*(xVec - theta));
    iter = 0;
    while (abs(f) > 1e-6)
        iter = iter + 1;
        if (f > 0)
            thetaMin = theta;
            theta = (theta + thetaMax)/2;
        else
            thetaMax = theta;
            theta = (thetaMin + theta)/2;
        end
        f = sum((q.^(exp(-theta)*theta.^xVec./xFactorial)).*(xVec - theta));
    end
    % f = sum((q.^(exp(-theta)*theta.^xVec./xFactorial)).*(xVec - theta));
end