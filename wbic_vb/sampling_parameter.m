function [parameter_ratio, parameter_mean] = sampling_parameter( phi, beta, u )
%sampling_vb �ϕ����㕪�z����p�����[�^���T���v�����O
%  phi:������̃n�C�p�[�p�����[�^(1*component), beta:�R���|�[�l���g�̕��U�̋t��(1*component),
%  u:�R���|�[�l���g�̒��S(component*dimension)
[component, dimension] = size(u);
PARAMETER_NUM = 10000;
parameter_ratio = zeros(PARAMETER_NUM, component);
parameter_mean = zeros(component,dimension, PARAMETER_NUM);
% disp(ones(dimension,1) * sqrt(1 ./ beta));
%%������,�R���|�[�l���g�̒��S��
for i = 1:1:PARAMETER_NUM
    parameter_ratio(i,:) = dirrnd(phi);
    parameter_mean(:,:,i) = sqrt(1 ./ beta)' * ones(1,dimension) .* randn(component,dimension) + u;
%     disp(parameter_mean(:,:,i));
end

end

