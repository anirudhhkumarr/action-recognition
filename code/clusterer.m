if exist('vl_version') ~= 3, run('vlfeat/toolbox/vl_setup') ; end

data = uint8(csvread('train-hoghof.csv',2));
K = 1000;
C = vl_ikmeans(data',K,'verbose');
C = C';
save('vocab_1k.mat', 'C');
