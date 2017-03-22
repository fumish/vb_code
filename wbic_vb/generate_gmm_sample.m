function x = generate_gmm_sample( ratio, mean, sample_number )
%generate_gmm_sample 混合正規分布に従うサンプルを作る
%   seed:乱数の種, ratio:真の混合比, mean:真の中心(component*dimension), sample_number:サンプル数

[component, dimension] = size(mean);

x = zeros(sample_number, dimension);
for i = 1:1:sample_number    
    index = rand;
    tmp = 0;
    for k = 1:1:component
        if tmp < index && index < tmp + ratio(k)
            x(i,:) = randn(1,dimension) + mean(k,:);
        end
        tmp = tmp + ratio(k);
    end
end

end

