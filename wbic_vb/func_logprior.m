function logprior = func_logprior(phi, beta, u, ratio, mean)
%%�n�C�p�[�p�����[�^��phi,beta,u�̎��O���z��
%%�p�����[�^�̒l��ratio,mean�̂Ƃ��̑ΐ��l���v�Z����
    dimension = size(u,2);
	logprior = gammaln(sum(phi))-sum(gammaln(phi))+sum((phi-1).*log(ratio))+ ...
	sum(dimension/2 * log(beta/(2*pi))) -sum(beta.*diag((mean-u)*(mean-u)')')/2;
end