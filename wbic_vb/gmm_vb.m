function [s,t,u, energy] = gmm_vb(samples, init_ratio, init_mean, beta, phi, inverse_temparature)
%%gmm_vb:変分ベイズを混合正規分布でやる
%%sample:サンプル(n*dimension) init_ratio:混合比の初期値(1*component)
%%init_mean:中心の初期値(component*dimension)
%%hyperparameter:ハイパーパラメータ([beta,phi])

[component, dimension] = size(init_mean);
n = size(samples,1);

ITERATION = 100;
s = n * init_ratio + phi;
t = n * init_ratio + beta;
u = init_mean;
% contribution = zeros(n,component);
for loop = 1:1:ITERATION
    %%Eステップ
    L = ones(n,1)*((psi(s)-psi(sum(s)))-diag(u*u')'/2-dimension./(2*t)) + samples*u' -diag(samples*samples')*ones(1,component)/2 + dimension/2*log(2*pi);
    L = inverse_temparature * L;
    maxL = max(L,[],2);
    L = L - maxL * ones(1,component);
    contribution = exp(L) ./ (sum(exp(L),2) * ones(1,component));    
    
    %%Mステップ
    s = inverse_temparature*sum(contribution,1) + phi;
    u = inverse_temparature*contribution' * samples ./ (t'*ones(1,dimension));
    t = inverse_temparature*sum(contribution,1) + phi;
end
F1 = sum((s-phi).*psi(s) - gammaln(s)) - n*psi(sum(s))+gammaln(sum(s))+sum(gammaln(phi))-gammaln(sum(phi));
F2 = sum(dimension/2*log(t./beta) + dimension*beta./(2*t) + diag(u*u')'.*beta/2)-component*dimension/2;
% F1 = sum((s-phi).*(psi(s)-psi(sum(s))) + dimension*beta./(2*t) + diag(u*u')'.*beta/2);
% F2 = sum(gammaln(phi))-gammaln(sum(phi))+dimension/2*sum(log(beta/(2*pi)));
% F3 = -sum(gammaln(s))+gammaln(sum(s))-dimension/2*sum(log(t/(2*pi)));
F3 = -sum(log(sum(exp(L),2)) + maxL);
energy = F1+F2+F3;

%%mixing_ratio = s / sum(s);
%%mean = u;

end