function ratio = dirrnd( phi )
%UNTITLED この関数の概要をここに記述
%   詳細説明をここに記述
    x = gamrnd(phi,ones(1,size(phi,2)));
    ratio = x ./ (sum(x));
end

