function [parameter_ratio, parameter_mean] = sampling_parameter( phi, beta, u )
%sampling_vb 変分事後分布からパラメータをサンプリング
%  phi:混合比のハイパーパラメータ(1*component), beta:コンポーネントの分散の逆数(1*component),
%  u:コンポーネントの中心(component*dimension)
[component, dimension] = size(u);
PARAMETER_NUM = 10000;
parameter_ratio = zeros(PARAMETER_NUM, component);
parameter_mean = zeros(component,dimension, PARAMETER_NUM);
% disp(ones(dimension,1) * sqrt(1 ./ beta));
%%混合比,コンポーネントの中心で
for i = 1:1:PARAMETER_NUM
    parameter_ratio(i,:) = dirrnd(phi);
    parameter_mean(:,:,i) = sqrt(1 ./ beta)' * ones(1,dimension) .* randn(component,dimension) + u;
%     disp(parameter_mean(:,:,i));
end

end

