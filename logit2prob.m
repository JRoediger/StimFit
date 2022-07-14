function probs = logit2prob(logs)
    odds = exp(logs);
    probs = odds./(1+odds);
end