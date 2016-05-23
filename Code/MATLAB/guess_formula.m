ratio = [];
for iter = 1:100
    iter
    m = 5;
    k = iter*2;
    
    % total = 0;
    % for i = 0:(m^k - 1)
    %     x = i;
    %     r = 0;
    %     a = ones(1, m);
    %     tmp = 0;
    %     for j = 1:k
    %         r = mod(x, m);
    %         a(r+1) = -a(r+1);
    %         if (a(r+1) == -1)
    %             tmp = tmp + 1;
    %         else
    %             tmp = tmp - 1;
    %         end
    %         x = (x-r)/m;
    %     end
    %     if (tmp == 0)
    %         total = total + 1;
    %     end
    % end
    %
    % [total, factorial(k)*(m/2)^(k/2)]
    % factorial(k)*(m/2)^(k/2)/total
    
    total = 0;
    for i = 0:m
        total = total + nchoosek(m,i)*(m-2*i)^k;
    end
    total = total/2^m;
    
    ratio(iter) = m^(0.9*k)/total*2^m;
    
end