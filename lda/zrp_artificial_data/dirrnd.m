function ratio = dirrnd( phi )
%UNTITLED ���̊֐��̊T�v�������ɋL�q
%   �ڍא����������ɋL�q
    x = gamrnd(phi,ones(1,size(phi,2)));
    ratio = x ./ (sum(x));
end

