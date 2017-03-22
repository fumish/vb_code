function x = generate_gmm_sample( ratio, mean, sample_number )
%generate_gmm_sample �������K���z�ɏ]���T���v�������
%   seed:�����̎�, ratio:�^�̍�����, mean:�^�̒��S(component*dimension), sample_number:�T���v����

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

