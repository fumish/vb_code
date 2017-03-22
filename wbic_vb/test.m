clf;
phi = [20,20];
for i = 1:1:30
    alpha = dirrnd(phi);
    plot(alpha(1),0,'r.');
    xlim([0,1]);
    ylim([0,1]);
    hold on;
end