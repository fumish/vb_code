function [x_data, t_data] = func_lda_data( phi, hop, DATA_SET, time_length, cell_num, car_number, outputname, seed_data )
%func_lda_data LDAのデータを作る
% phi:コンポーネントの比率を決めるパラメータ(Dir(|phi)に従ってratioを発生させる)
% hop:各コンポーネントの前進確率(TRUE_K*M)
% time_length: 各データセットでのセルオートマトンの時間(1*DATA_SET)
% cell_num:各データセットでのセルの個数(1*DATA_SET)
% car_number:各データセットでの車の数(1*DATA_SET)
% seed:乱数の種
% outputname:結果を出力するファイル名

rng(seed_data);

[TRUE_K,M] = size(hop);
x_data = zeros(max(car_number),M,DATA_SET); %%car_numberはデータセットごとに可変にできるように最終的にがんばる
t_data = zeros(max(car_number),M,DATA_SET);
% test_d_ratio = dirrnd(phi);
ratio = zeros(TRUE_K,DATA_SET);
for d = 1:1:DATA_SET
    %%各コンポーネントの発生確率を生成
    d_ratio = dirrnd(phi);
%     d_ratio = test_d_ratio;
    
    ratio(:,d) = d_ratio';
    [d_x_data, d_t_data, d_label, ~] = func_zrp_simulation_main_ring(d_ratio, hop, M, car_number(d), cell_num(d), time_length(d));
    x_data(1:car_number(d),:,d) = d_x_data;
    t_data(1:car_number(d),:,d) = d_t_data;
    
    label_prop = zeros(1,TRUE_K);
    for k = 1:1:TRUE_K
        arg = find(d_label == k);
        label_prop(k) = length(arg)/length(d_label);
    end
    
%     dlmwrite(strcat('true_', outputname), d_ratio', '-append');
    if d == 1
        dlmwrite(strcat('feature_', outputname), [d_x_data, d_t_data, d_label, d_x_data ./ d_t_data]);
    else
        dlmwrite(strcat('feature_', outputname), [d_x_data, d_t_data, d_label, d_x_data ./ d_t_data], '-append', 'roffset', 1);
    end
    dlmwrite(strcat('feature_', outputname), label_prop, '-append');
end
% disp(x_data);
% disp(t_data);
dlmwrite(strcat('true_', outputname), ratio);
dlmwrite(strcat('true_', outputname), hop, '-append', 'roffset', 1);

end

