function fundamental_diagram_main( phi, hop, time_length, time_division, cell_num, seed, outputname )
%%% fundamental_diagram_main LDA-TASEP�Ɋ�Â����V�~�����[�V�����ɂ���{�}���v�Z
%%% �e���������Ƃɗ��ʂ��o��
%%% --------------�����Ɩ߂�l-------------------
%%% phi:�g�s�b�N�����𐧌䂷��f�B���N�����z�̃p�����[�^(1*K)
%%% hop:�O�i�m��(1*K)
%%% time_length:�S����
%%% time_division:���Ԃ̕�����
%%% cell_num:�Z���̐�
%%% seed:�����̎�
%%% ---------------------------------------------

% rng(seed);
fd_dir = './fundamental_diagram/';
mkdir(fd_dir);

feature_dir = './feature/';
mkdir(feature_dir);

common_file = strcat('cell', num2str(cell_num), '_time', num2str(time_length), '_', outputname);

INIT_NUM = 100; %%�����ԗ����x�ɑ΂��ĕς��鏉���l�̉�
% INIT_NUM = 1;
rho = 0.05:0.05:0.95; %%�Z���̌��ɑ΂���ԗ����x

for i = 1:1:length(rho)
    %%�e���x�ɑ΂��ė��ʂ̌v�Z���s��
    car_number = floor(cell_num * rho(i));
    
    for l = 1:1:INIT_NUM
        rng(seed+(i-1)*INIT_NUM+l);
        %%�f�[�^�������s��
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
            
            %%���ʂ̏o��
            if i == 1 && l == 1
                dlmwrite(strcat(fd_dir, 'fd_', num2str(d), 'division_' , common_file), [rho(i), d_flux]);
            else
                dlmwrite(strcat(fd_dir, 'fd_', num2str(d), 'division_' , common_file), [rho(i), d_flux], '-append');
            end
        end
        
        %%�����ʂ̏o��
        dlmwrite(strcat(feature_dir, 'feature_', num2str(car_number), 'car_', num2str(seed+(i-1)*INIT_NUM+l), 'seed_', common_file), x_data );
        dlmwrite(strcat(feature_dir, 'feature_', num2str(car_number), 'car_', num2str(seed+(i-1)*INIT_NUM+l), 'seed_', common_file), t_data, '-append');
    end
end

end

