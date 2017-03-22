function [x_data, t_data, label, final_flux, final_car_pos] = func_tasep_simulation_main_ring(init_car_pos, ratio, hop, car_num, cell_num,  time_length)
%%%tomei_simulation_main: ������𕡐��̃N���X�̎ԗ��������V�~�����[�V���������
%%%po:�����m�� raito:������(1*K), hop:�O�i�m��(K*3), cell_num:�Z���̐�,
%%%total_time:�V�~�����[�V�����̑S����, init,last_cell:���ڃZ���̑O��

COMPONENT = size(ratio,2);
%%�Ԃ̏����擾

%%�����ݒ�
car_information = zeros(car_num,3);
x_data = zeros(car_num,1);
t_data = zeros(car_num,1);
%����(i,:)��i�Ԗڂ̎Ԃ̏��
%����1�ڂ����݂̎Ԃ̈ʒu(�o����-1�ɂ���),
%2�`vec_size+1�ڂ��O�i�m��,vec_size+2�ڂ��N���X�ԍ�
cell_automaton = zeros(time_length,cell_num);
%�����Ɋe�������Ƃ̃Z���I�[�g�}�g���̏�Ԃ�����
%0�͉������Ȃ����Ƃ�\���Ă���

%%%�Z���I�[�g�}�g���̏�����Ԃ�ݒ�
for i = 1:1:car_num
    cell_automaton(1,floor(i*cell_num/car_num)) = i;
    
    tmp_r = rand;
    ind = 0;
    for k = 1:1:COMPONENT
        if ind < tmp_r && tmp_r < ind + ratio(k)
            %%k�Ԗڂ̃z�b�v�m�������ԗ��𔭐�������
            if isempty(init_car_pos)
                car_information(i, 1) = floor(i*cell_num/car_num);
            else
                car_information(i,1) = init_car_pos(i);
            end
            car_information(i, 2) = hop(k);
            car_information(i, 3) = k;
            break;
        end
        ind = ind + ratio(k);
    end
end


%%�V�~�����[�V�����J�n
flux = zeros(time_length-1,1);
for t = 1:1:time_length-1
    %%�Z���I�[�g�}�g���̍X�V
    current_pos = car_information(:,1);
    for i = 1:1:car_num
        %%�O�̎ԂƂ̎Ԋԋ���
        if i == car_num
            delta_x = mod((car_information(1,1)-car_information(car_num,1))-1,cell_num);
        else
            delta_x = mod((car_information(i+1,1)-car_information(i,1))-1,cell_num);
        end
        if delta_x == 0
            cell_automaton(t+1,current_pos(i)) = i;
            continue;
        else
            if rand < car_information(i,2)
                %%�O�ɐi��
                flux(t) = flux(t) + 1;
                x_data(i,1) = x_data(i,1) + 1;
                t_data(i,1) = t_data(i,1) + 1;
                if current_pos(i) == cell_num
                    cell_automaton(t+1,1) = i;
                    car_information(i,1) = 1;
                else
                    cell_automaton(t+1, current_pos(i)+1) = i;
                    car_information(i,1) = current_pos(i)+1;
                end
            else
                cell_automaton(t+1, current_pos(i)) = i;
                t_data(i,1) = t_data(i,1) + 1;
            end            
        end
    end    
end
label = car_information(:,3);
final_flux = flux(length(flux)) / cell_num;
final_car_pos = car_information(:,1);
end