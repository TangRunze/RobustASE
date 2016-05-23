function errorRate = errorratecalculator(tauStar, tauHat, nBlock)

nv = (tauStar == 0);
tauStar(nv) = [];
tauHat(nv) = [];
n = length(tauStar);

errorRate = n;

permutation = perms(1:nBlock);
for iFactorial = 1:factorial(nBlock)
    pos = permutation(iFactorial,:);
    tauTmp = tauHat;
    for jBlock = 1:nBlock
        nv = (tauHat == pos(jBlock));
        tauTmp(nv) = jBlock;
    end
    if sum(tauStar ~= tauTmp) < errorRate
        errorRate = sum(tauStar ~= tauTmp);
    end
end

errorRate = errorRate/n;