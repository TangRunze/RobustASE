function theta = lqexpectation_exp(theta0, q, c, epsilon)
    aTmp = (-c*epsilon*(-1+q)^2*q+(-1+epsilon)*(-1+q)^2*q*theta0);
    bTmp = (c^2*epsilon*(-1+q)^2+2*c*(-1+q)*q*theta0+(1-epsilon)*(-1+q)^2*theta0^2);
    cTmp = (-c^2*(epsilon*(-2+q)+q)*theta0-c*(-2+2*epsilon+2*q-epsilon*q)*theta0^2);
    dTmp = c^2*theta0^2;
    pTmp = -bTmp/3/aTmp;
    qTmp = pTmp^3+(bTmp*cTmp-3*aTmp*dTmp)/(6*aTmp^2);
    rTmp=cTmp/3/aTmp;
    theta = nthroot(qTmp+sqrt(qTmp^2+(rTmp-pTmp^2)^3), 3) ...
        + nthroot(qTmp-sqrt(qTmp^2+(rTmp-pTmp^2)^3), 3) ...
        + pTmp;
end