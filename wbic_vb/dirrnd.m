function alpha = dirrnd( phi )
%dirrnd �f�B���N�����z����̃T���v���𐶐�
%   �ڍא����������ɋL�q

x = gamrnd(phi, ones(1,size(phi,2)));
alpha = x / sum(x);

end

