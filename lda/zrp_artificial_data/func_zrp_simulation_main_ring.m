function [x_data, t_data, label, final_flux] = func_zrp_simulation_main_ring(ratio, hop, vec_size, car_num, cell_num,  time_length)
%%%tomei_simulation_main: 東名上を複数のクラスの車両が動くシミュレーションを作る
%%%po:発生確率 raito:混合比(1*K), hop:前進確率(K*3), cell_num:セルの数,
%%%total_time:シミュレーションの全時間, init,last_cell:注目セルの前後

COMPONENT = size(ratio,2);
%%車の情報を取得

%%初期設定
car_information = zeros(car_num,vec_size+2);
x_data = zeros(car_num,vec_size);
t_data = zeros(car_num,vec_size);
%この(i,:)がi番目の車の情報
%情報は1つ目が現在の車の位置(出たら-1にする),
%2～vec_size+1つ目が前進確率,vec_size+2つ目がクラス番号
cell_automaton = zeros(time_length,cell_num);
%ここに各時刻ごとのセルオートマトンの状態を入れる
%0は何もいないことを表している

%%%セルオートマトンの初期状態を設定
for i = 1:1:car_num
    cell_automaton(1,floor(i*cell_num/car_num)) = i;
    
    tmp_r = rand;
    ind = 0;
    for k = 1:1:COMPONENT
        if ind < tmp_r && tmp_r < ind + ratio(k)
            %%k番目のホップ確率を持つ車両を発生させる
            car_information(i, 1) = floor(i*cell_num/car_num);
            car_information(i, 2:vec_size+1) = hop(k,:);
            car_information(i, vec_size+2) = k;
            break;
        end
        ind = ind + ratio(k);
    end
end


%%シミュレーション開始
flux = zeros(time_length-1,1);
for t = 1:1:time_length-1
    %%セルオートマトンの更新
    current_pos = car_information(:,1);
    for i = 1:1:car_num
        %%前の車との車間距離
        if i == car_num
            delta_x = mod((car_information(1,1)-car_information(car_num,1))-1,cell_num);
        else
            delta_x = mod((car_information(i+1,1)-car_information(i,1))-1,cell_num);
        end
        if delta_x > vec_size
            delta_x = vec_size;
        end
        if delta_x == 0
            cell_automaton(t+1,current_pos(i)) = i;
            continue;
        else
            if rand < car_information(i,1+delta_x)
                %%前に進む
                flux(t) = flux(t) + 1;
                x_data(i,delta_x) = x_data(i,delta_x) + 1;
                t_data(i,delta_x) = t_data(i,delta_x) + 1;
                if current_pos(i) == cell_num
                    cell_automaton(t+1,1) = i;
                    car_information(i,1) = 1;
                else
                    cell_automaton(t+1, current_pos(i)+1) = i;
                    car_information(i,1) = current_pos(i)+1;
                end
            else
                cell_automaton(t+1, current_pos(i)) = i;
                t_data(i,delta_x) = t_data(i,delta_x) + 1;
            end            
        end
    end    
end
label = car_information(:,vec_size+2);
final_flux = flux(time_length-1) / cell_num;
end