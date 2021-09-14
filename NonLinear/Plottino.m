% Plottino.m

figure(1);
plot(x_out1.data.*cos(psi.data),y_out1.data.*sin(psi.data))
grid on;
xlabel("p_x");
ylabel("p_y");