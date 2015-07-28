function theta = l1expectation_exp(theta0, c, epsilon)
    theta = - c*theta0/(c*epsilon - c - epsilon*theta0);
end