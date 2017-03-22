function [min_energy, min_ratio, min_hop, min_contribution] = func_vb_lda_tasep(x_data, t_data, car_number, COMPONENT, phi_val, alpha_val, seed_learn)
%%%%%EM: tasep_dataにVBアルゴリズムをやる
%%%%% seed:乱数の種 parameter:mixing_ratio(COMPONENT * DATA_SET),
%%%%% hop_prob(COMPONENT*v_len),
%%%%% contribution:寄与率(car_num * COMPONENT * DATA_SET) v_len:考慮する車間距離

rng(seed_learn);

ITERATION = 1000;
[N, DATA_SET] = size(x_data);
%n = size(zrp_data,1);

% inverse_temparature = 1;
phi = phi_val * ones(COMPONENT, DATA_SET);
alpha = alpha_val * ones(COMPONENT, 1);
beta = alpha_val * ones(COMPONENT, 1);

min_energy = realmax;
SEED_NUM = 100;
for l = 1:1:SEED_NUM
    
    contribution = zeros(N, COMPONENT, DATA_SET);
    %%LDAの隠れ変数の初期化
    for d = 1:1:DATA_SET
        for n = 1:1:car_number(d)
            contribution(n,:,d) = dirrnd(ones(1,COMPONENT));
        end
    end
    
    %%変分ベイズをやるとこ
    for t = 1:1:ITERATION
        %%Mステップ  
        s = squeeze(sum(contribution,1)) + phi;

        T = zeros(COMPONENT,1);
        U = zeros(COMPONENT,1);
        for d = 1:1:DATA_SET
            T = T + contribution(1:car_number(d),:,d)'*x_data(1:car_number(d),d) + alpha;
            U = U + contribution(1:car_number(d),:,d)'*(t_data(1:car_number(d),d)-x_data(1:car_number(d),d)) + beta;
        end
        
        %%Eステップ
        L = zeros(size(contribution));
        for d = 1:1:DATA_SET           
            L(1:car_number(d),:,d) = ones(car_number(d),1)*psi(s(:,d))'+ones(car_number(d),COMPONENT)*psi(sum(s(:,d))) ...
                +x_data(1:car_number(d),d)*psi(T)' ...
                +(t_data(1:car_number(d),d)-x_data(1:car_number(d),d))*psi(U)' - t_data(1:car_number(d),d)*psi(T+U)' ...
                +(gammaln(t_data(1:car_number(d),d)+1)-gammaln(x_data(1:car_number(d),d)+1) ...
                -gammaln(t_data(1:car_number(d),d)-x_data(1:car_number(d),d)+1))*ones(1,COMPONENT);
        end
        for d = 1:1:DATA_SET
            d_maxL = max(L(1:car_number(d),:,d),[],2);
            contribution(1:car_number(d),:,d) = exp(L(1:car_number(d),:,d)-d_maxL * ones(1,COMPONENT)) ./ (sum(exp(L(1:car_number(d),:,d) - d_maxL * ones(1,COMPONENT)),2) * ones(1,COMPONENT));
        end
        
    end
    %%エネルギー計算
    F1 = sum(sum( (s - phi) .* ( psi(s) - ones(COMPONENT,1)*psi(sum(s,1)) ) ) );
    F2 = sum( (T-alpha) .* (psi(T)-psi(T+U)) + (U-alpha) .* (psi(U)-psi(T+U)));
    F3 = -sum(sum(gammaln(s),1) - gammaln(sum(s,1))) - sum(gammaln(T) + gammaln(U) - gammaln(T+U));
    F4 = sum(sum(gammaln(phi),1) - gammaln(sum(phi,1))) + sum(gammaln(alpha) + gammaln(alpha) - gammaln(alpha+alpha));
    F5 = 0;
    for d = 1:1:DATA_SET
        d_maxL = max(L(1:car_number(d),:,d),[],2);
        F5 = F5 -sum(log(sum(exp(L(1:car_number(d),:,d) - d_maxL * ones(1,COMPONENT)),2)) + d_maxL);
%         disp([L(:,:,d), contribution(:,:,d)]);
    end
    
    F = F1 + F2 + F3 + F4 + F5;
%     disp([F1+F2+F3+F4,F5]);
    if F < min_energy
        min_energy = F;
        min_contribution = contribution;
        min_ratio = s ./ (ones(COMPONENT,1)*sum(s,1));
        min_hop = T ./ (T + U);
    end
%     disp(l);
%     disp(F);
end

end







