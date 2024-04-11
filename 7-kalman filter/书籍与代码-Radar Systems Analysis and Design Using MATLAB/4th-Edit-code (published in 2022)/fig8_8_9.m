
figure (5)
plot(doppler,x(:,149),'b', 'linewidth',1);
xlabel ('\bf Frequency - Hz')
ylabel ('\bf Ambiguity - Vots')
grid on
title('\bf \chi (0, f) - Zero delay cut')
figure (6)
plot(time,x(228,:),'b','linewidth',1);
xlabel ('\bf delay - seconds')
ylabel ('\bf Ambiguity - Vots')
grid on
title('\bf \chi (\tau,0) - Zero Doppler cut')