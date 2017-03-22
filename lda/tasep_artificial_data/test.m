clear;
phi = [10,10,10];
hop = [0.2,0.5,0.7];
% hop = [0.8,0.9,0.9; 0.4,0.7,0.7; 0.1,0.7,0.9; 0.1,0.3,0.7];
DATA_SET = 6;

time_length = 500*ones(1,DATA_SET);
cell_num = 200*ones(1,DATA_SET);
% car_number = 100*ones(1,DATA_SET);
car_number = [50,60,70,80,90,100];

seed_data = 1;
seed_learn = 2;
outputname = '0826_test.csv';

LEARNING_COMPONENT = 4;
phi_const = 0.1;
alpha_const = 1;

vb_lda_main(phi, hop, DATA_SET, time_length, cell_num, car_number, seed_data, seed_learn, LEARNING_COMPONENT, phi_const, alpha_const, outputname);
% func_lda_data(phi, hop, DATA_SET, time_length, cell_num, car_number, seed_data, outputname);