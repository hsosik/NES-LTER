
function pVal = phaseScrambling(Y)
    % Adapted from Solow and Duplisea (2007)
    % Y is the entire timeseries of each component j (i.e., species, spatial location, etc.)
    % m is the number of components j=1,..m
    [n,m] = size(Y);
    R = variance_ratio(Y);
    Y_star = zeros(n,m);
    rep = 1000;
    R_star = zeros(rep,1);

    for r=1:rep
        for j=1:m
            Yj_bar = nanmean(Y(:,j)); 
            
            for t=1:n
                summation=0;
                for k=1:n
                    for l=1:n
                        Ul = 2*pi*rand(1,1);
                        summation = summation + (Y(k,j) - Yj_bar)*cos(2*pi*(l-1)*(k-t)/n + Ul);
                    end
                end
                Y_star(t,j) = Yj_bar + (sqrt(2)/n)*summation;
            end
        end

        R_star(r) = variance_ratio(Y_star);
    end
    
    % Testing the null hypothesis of independence, i.e., R = 1
    if R < 1
        pVal = sum(R_star < R)/rep;
    else
        pVal = sum(R_star > R)/rep;
    end
end