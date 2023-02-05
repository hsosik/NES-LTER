function vr = variance_ratio(input_data)
    % input_data matrix has taxa columns, sample rows
    % vr output is variance ratio from Schluter 19
    % following full calc from Schluter, Eq 2,3,5
    % Heidi M. Sosik, Woods Hole Oceanographic Institution, Jan 2023
    
    %eliminate any rows with NaNs
    input_data(isnan(sum(input_data,2)),:) = [];
    
    if size(input_data,1) > 1
            T = sum(input_data,2);
        t = mean(T);
        s_T2_N = sum((T-t).^2); %1/N from paper cancels in final ratio
        t_i = mean(input_data);
        t_imat = repmat(t_i,size(input_data,1),1);
        sigma_i2_N = sum((input_data-t_imat).^2);
        vr(1) = s_T2_N/sum(sigma_i2_N);
    else
        vr(1) = NaN; %case for only 1 observation
    end    
    
    %matlab version of codyn
    %This works too BUT it's 3 times slower than above!!
   % vr(1) = sum(sum(cov(input_data))/sum(var(input_data)),2);

    %Isabel Honda's approach
    %%SEEMS Wrong or is it a different Variance Ratio than Schluter??
%    s_T2 = sum(var(input_data),2); 
%    covSum = sum(cov(input_data),'all');
%    varTest = (s_T2+2*covSum)/s_T2;
%    vr(2) = varTest;
 
%     %Ji's approach
%     varsum = sum(var(input_data));
%     covsum = sum(triu(cov(input_data),1),'all');
%     vr(3)=(varsum+2*covsum)/varsum;
%     
%     %codyn approach
%     vr(4) = sum(cov(input_data), 'all')/sum(var(input_data, 'omitnan'),2);
    
end

