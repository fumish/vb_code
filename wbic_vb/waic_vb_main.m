function wbic_vb = waic_vb_main(true_ratio, true_mean, sample_number, LEARNING_COMPONENT, phi_const, beta_const, seed)
%%変分事後分布を提案分布としたwbicの計算
%%trueとlearnerが混合正規分布でやってみる
%%true_ratio:真の混合比 true_mean:コンポーネントの真の中心
%%sample_number:サンプル数, LEARNING_COMPONENT:learnerのコンポーネント数
%%seed:乱数の種 phi,beta:ハイパーパラメータ

%%乱数初期化
rng(seed);

dimension = size(true_mean, 2);
% inverse_temparature = 1 / log(sample_number);
inverse_temparature = 1;

phi = phi_const * ones(1,LEARNING_COMPONENT);
beta = beta_const * ones(1, LEARNING_COMPONENT);

%%サンプルを発生させる
samples = generate_gmm_sample(true_ratio, true_mean, sample_number);
mean = sum(samples,1) / sample_number;
%%変分ベイズで学習を行う
VB_TIME = 100;
min_energy = realmax;
min_phi = zeros(1,LEARNING_COMPONENT);
min_mean = zeros(LEARNING_COMPONENT, dimension);
min_v = zeros(1,LEARNING_COMPONENT);
for t = 1:1:VB_TIME
	%%初期値を決める
	init_ratio = dirrnd(ones(1,LEARNING_COMPONENT));
	init_mean = ones(LEARNING_COMPONENT,1)*mean + 0.01*randn(LEARNING_COMPONENT, dimension);
	
	[est_phi, est_v, est_mean, est_energy] = gmm_vb(samples, init_ratio, init_mean, beta, phi, inverse_temparature);
	if min_energy > est_energy
		min_energy = est_energy;
		min_phi = est_phi;
		min_mean = est_mean;
		min_v = est_v;
%         disp(min_phi / sum(min_phi));
%         disp(min_mean);
	end
end

%%変分ベイズからサンプルを取得
% disp(min_energy);
% disp(min_energy + func_loglikelihood(samples,true_ratio,true_mean)); 
disp(min_phi / sum(min_phi));
% disp(min_phi);
disp(min_mean);
[parameter_ratio, parameter_mean] = sampling_parameter(min_phi, min_v, min_mean);
plot(parameter_ratio(:,1),parameter_ratio(:,2),'r.');
% disp(parameter_ratio);
% disp(parameter_mean);
% wbic_vb = 0;
%%提案分布と事後分布の比を計算
PARAMETER_NUM = size(parameter_ratio,1);
%sampling_ratio = zeros(PARAMETER_NUM,1);
bayes_value = zeros(PARAMETER_NUM,1);
vb_value = zeros(PARAMETER_NUM,1);
for i = 1:1:PARAMETER_NUM
    bayes_value(i) = inverse_temparature * func_loglikelihood(samples,parameter_ratio(i,:),parameter_mean(:,:,i)) + func_logprior(phi,beta,zeros(LEARNING_COMPONENT,dimension),parameter_ratio(i,:),parameter_mean(:,:,i));
    vb_value(i) = func_logprior(min_phi, min_v, min_mean, parameter_ratio(i,:), parameter_mean(:,:,i));
end
sampling_ratio = exp(bayes_value-vb_value);
format long;
dlmwrite('sampling_ratio.csv',sampling_ratio);

%%wbic_vbをやる
partial_wbic_vb = zeros(PARAMETER_NUM,1);
for i = 1:1:PARAMETER_NUM
	partial_wbic_vb(i) = -func_loglikelihood(samples,parameter_ratio(i,:),parameter_mean(:,:,i)) * sampling_ratio(i);
end
wbic_vb = sum(partial_wbic_vb) / sum(sampling_ratio) + func_loglikelihood(samples, true_ratio, true_mean);

end