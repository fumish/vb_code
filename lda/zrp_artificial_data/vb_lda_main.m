function vb_lda_main( phi, hop, DATA_SET, time_length, cell_num, car_number, seed_data, seed_learn, LEARNING_COMPONENT, phi_const, alpha_const, outputname)
%vb_zrp_main 周期状のセルオートマトンの推定
%  

%%データ発生
[x_data, t_data] = func_lda_data(phi, hop, DATA_SET, time_length, cell_num, car_number, outputname, seed_data);
[min_energy, min_ratio, min_hop, min_contribution] = func_vb_lda(x_data, t_data, LEARNING_COMPONENT, phi_const, alpha_const, seed_learn);


for d = 1:1:DATA_SET
    if d == 1
        dlmwrite(strcat('contribution_', outputname), min_contribution(:,:,d));
    else
        dlmwrite(strcat('contribution_', outputname), min_contribution(:,:,d), '-append', 'roffset', 1);
    end
end
dlmwrite(strcat('estparameter_', outputname), min_ratio);
dlmwrite(strcat('estparameter_', outputname), min_hop, '-append', 'roffset', 1);

disp(min_energy);

end

