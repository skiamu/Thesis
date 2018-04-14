function [VaR,ES] = hHistoricalVaRES(Sample,VaRLevel)
    % Compute historical VaR and ES
    % See [4] for technical details

    % Convert to losses
    Sample = -Sample;
    
    N = length(Sample);
    k = ceil(N*VaRLevel);
    
    z = sort(Sample);
    
    VaR = z(k);
    
    if k < N
       ES = ((k - N*VaRLevel)*z(k) + sum(z(k+1:N)))/(N*(1 - VaRLevel));
    else
       ES = z(k);
    end
end

function [VaR,ES] = hNormalVaRES(Mu,Sigma,VaRLevel)
    % Compute VaR and ES for normal distribution
    % See [3] for technical details
    
    VaR = -1*(Mu-Sigma*norminv(VaRLevel));
    ES = -1*(Mu-Sigma*normpdf(norminv(VaRLevel))./(1-VaRLevel));

end

function [VaR,ES] = hTVaRES(DoF,Mu,Sigma,VaRLevel)
    % Compute VaR and ES for t location-scale distribution
    % See [3] for technical details

    VaR = -1*(Mu-Sigma*tinv(VaRLevel,DoF));
    ES_StandardT = (tpdf(tinv(VaRLevel,DoF),DoF).*(DoF+tinv(VaRLevel,DoF).^2)./((1-VaRLevel).*(DoF-1)));
    ES = -1*(Mu-Sigma*ES_StandardT);

end