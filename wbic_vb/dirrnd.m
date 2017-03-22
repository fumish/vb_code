function alpha = dirrnd( phi )
%dirrnd ディリクレ分布からのサンプルを生成
%   詳細説明をここに記述

x = gamrnd(phi, ones(1,size(phi,2)));
alpha = x / sum(x);

end

