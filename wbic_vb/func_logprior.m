function logprior = func_logprior(phi, beta, u, ratio, mean)
%%ハイパーパラメータがphi,beta,uの事前分布の
%%パラメータの値がratio,meanのときの対数値を計算する
    dimension = size(u,2);
	logprior = gammaln(sum(phi))-sum(gammaln(phi))+sum((phi-1).*log(ratio))+ ...
	sum(dimension/2 * log(beta/(2*pi))) -sum(beta.*diag((mean-u)*(mean-u)')')/2;
end