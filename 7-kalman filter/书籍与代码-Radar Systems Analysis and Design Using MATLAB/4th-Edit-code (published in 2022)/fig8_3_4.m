
figure (5)
plot(fdy,x(:,76),'b', 'linewidth',1);
xlabel ('\bf Frequency - Hz')
ylabel ('\bf Ambiguity - Vots')
grid on
title('\bf \chi (0, f) - Zero delay cut')
figure (6)
plot(taux,x(76,:),'b','linewidth',1);
xlabel ('\bf delay - seconds')
ylabel ('\bf Ambiguity - Vots')
grid on
title('\bf \chi (\tau,0) - Zero Doppler cut')