function [min_energy, min_ratio, min_hop, min_contribution] = func_vb_lda(x_data, t_data, COMPONENT, phi_val, alpha_val, seed_learn)
%%%%%EM: tasep_data��VB�A���S���Y�������
%%%%% seed:�����̎� parameter:mixing_ratio(COMPONENT * DATA_SET),
%%%%% hop_prob(COMPONENT*v_len),
%%%%% contribution:��^��(car_num * COMPONENT * DATA_SET) v_len:�l������Ԋԋ���

rng(seed_learn);

ITERATION = 1000;
[N, M, DATA_SET] = size(x_data);
%n = size(zrp_data,1);

% inverse_temparature = 1;
phi = phi_val * ones(COMPONENT, DATA_SET);
alpha = alpha_val * ones(COMPONENT, M);
beta = alpha_val * ones(COMPONENT, M);

min_energy = realmax;
SEED_NUM = 100;
for l = 1:1:SEED_NUM
    
    contribution = zeros(N, COMPONENT, DATA_SET);
    %%LDA�̉B��ϐ��̏�����
    for n = 1:1:N
        for d = 1:1:DATA_SET
            contribution(n,:,d) = dirrnd(ones(1,COMPONENT));
        end
    end
    
    %%�ϕ��x�C�Y�����Ƃ�
    for t = 1:1:ITERATION
        %%M�X�e�b�v
        s = squeeze(sum(contribution,1)) + phi;
        T = zeros(COMPONENT,M);
        U = zeros(COMPONENT,M);
        for d = 1:1:DATA_SET
            T = T + contribution(:,:,d)'*x_data(:,:,d) + alpha;
            U = U + contribution(:,:,d)'*(t_data(:,:,d)-x_data(:,:,d)) + beta;
        end
        
        %%E�X�e�b�v
        L = zeros(size(contribution));
        for d = 1:1:DATA_SET
            L(:,:,d) = sum((gammaln(t_data(:,:,d)+1)-gammaln(x_data(:,:,d)+1)-gammaln(t_data(:,:,d)-x_data(:,:,d)+1)),2)*ones(1,COMPONENT) + ones(N,1)*psi(s(:,d))'+x_data(:,:,d)*psi(T)'+(t_data(:,:,d)-x_data(:,:,d))*psi(U)'-t_data(:,:,d)*psi(T+U)';
        end
        for d = 1:1:DATA_SET
            d_maxL = max(L(:,:,d),[],2);
            contribution(:,:,d) = exp(L(:,:,d)-d_maxL * ones(1,COMPONENT)) ./ (sum(exp(L(:,:,d) - d_maxL * ones(1,COMPONENT)),2) * ones(1,COMPONENT));
        end
        
    end
    %%�G�l���M�[�v�Z
    F1 = sum(sum( (s - phi) .* (psi(s)-psi(ones(COMPONENT,1)*sum(s,1))) ));
    F2 = sum(sum( (T-alpha) .* (psi(T)-psi(T+U)) + (U-alpha) .* (psi(U)-psi(T+U)) ));
    F3 = -sum(sum(gammaln(s),1) - gammaln(sum(s,1))) - sum(sum(gammaln(T) + gammaln(U) - gammaln(T+U)));
    F4 = sum(sum(gammaln(phi),1) - gammaln(sum(phi,1))) + sum(sum(gammaln(alpha) + gammaln(alpha) - gammaln(alpha+alpha)));
    F5 = 0;
    for d = 1:1:DATA_SET
        d_maxL = max(L(:,:,d),[],2);
        F5 = F5 -sum(log(sum(exp(L(:,:,d) - d_maxL * ones(1,COMPONENT)),2)) + d_maxL);
        disp([L(:,:,d), contribution(:,:,d)]);
    end
    % disp(F1+F2+F3+F4);
    F = F1 + F2 + F3 + F4 + F5;

    if F < min_energy
        min_energy = F;
        min_contribution = contribution;
        min_ratio = s ./ (ones(COMPONENT,1)*sum(s,1));
        min_hop = T ./ (T + U);
    end
    disp(l);
    disp(F);
end

end







