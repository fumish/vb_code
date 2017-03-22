function loglikelihood = func_loglikelihood(samples, ratio, mean)
%%%�p�����[�^(ratio,mean)�̑ΐ��ޓx���v�Z����

n = size(samples,1);
[component, dimension] = size(mean);
partial = zeros(n,1);
for i = 1:1:n
comp_dist = diag((ones(component,1)*samples(i,:)-mean)*(ones(component,1)*samples(i,:)-mean)')';
partial(i,:) = log(sum(ratio .* exp(-comp_dist/2) / (sqrt(2*pi)^dimension)));
end

loglikelihood = sum(partial);
end
