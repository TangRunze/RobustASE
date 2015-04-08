function theta = lqsolve(xVec, q)
    theta = mean(xVec);
    thetaMin = min(xVec);
    thetaMax = max(xVec);
    xFactorial = factorial(xVec);
    f = sum((theta.^xVec./xFactorial).^(1 - q).*(xVec - theta));
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
        f = sum((theta.^xVec./xFactorial).^(1 - q).*(xVec - theta));
    end
    iter
end