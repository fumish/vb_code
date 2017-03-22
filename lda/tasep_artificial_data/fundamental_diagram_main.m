function fundamental_diagram_main( phi, hop, time_length, time_division, cell_num, seed, outputname )
%%% fundamental_diagram_main LDA-TASEPに基づいたシミュレーションによる基本図を計算
%%% 各時分割ごとに流量を出力
%%% --------------引数と戻り値-------------------
%%% phi:トピック生成を制御するディリクレ分布のパラメータ(1*K)
%%% hop:前進確率(1*K)
%%% time_length:全時間
%%% time_division:時間の分割数
%%% cell_num:セルの数
%%% seed:乱数の種
%%% ---------------------------------------------

% rng(seed);
fd_dir = './fundamental_diagram/';
mkdir(fd_dir);

feature_dir = './feature/';
mkdir(feature_dir);

common_file = strcat('cell', num2str(cell_num), '_time', num2str(time_length), '_', outputname);

INIT_NUM = 100; %%同じ車両密度に対して変える初期値の回数
% INIT_NUM = 1;
rho = 0.05:0.05:0.95; %%セルの個数に対する車両密度

for i = 1:1:length(rho)
    %%各密度に対して流量の計算を行う
    car_number = floor(cell_num * rho(i));
    
    for l = 1:1:INIT_NUM
        rng(seed+(i-1)*INIT_NUM+l);
        %%データ分割を行う
        current_car_pos = [];
        
        x_data = zeros(car_number, time_division);
        t_data = zeros(size(x_data));
        
        for d = 1:1:time_division
            d_ratio = dirrnd(phi);
            
            init_time = 1 + (d-1)*floor(time_length/cell_num);
            last_time = d*floor(time_length/cell_num);
            
            [d_x_data, d_t_data, ~, d_flux, current_car_pos] = func_tasep_simulation_main_ring(current_car_pos, d_ratio, hop, car_number, cell_num, last_time-init_time+1);
            x_data(:,d) = d_x_data;
            t_data(:,d) = d_t_data;            
            
            %%流量の出力
            if i == 1 && l == 1
                dlmwrite(strcat(fd_dir, 'fd_', num2str(d), 'division_' , common_file), [rho(i), d_flux]);
            else
                dlmwrite(strcat(fd_dir, 'fd_', num2str(d), 'division_' , common_file), [rho(i), d_flux], '-append');
            end
        end
        
        %%特徴量の出力
        dlmwrite(strcat(feature_dir, 'feature_', num2str(car_number), 'car_', num2str(seed+(i-1)*INIT_NUM+l), 'seed_', common_file), x_data );
        dlmwrite(strcat(feature_dir, 'feature_', num2str(car_number), 'car_', num2str(seed+(i-1)*INIT_NUM+l), 'seed_', common_file), t_data, '-append');
    end
end

end

