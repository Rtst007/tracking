function plot_figures_chap8( ambig, delay, freq)
% This function is used to plot figures in Chapter 6
%
figure
mesh(delay,freq,(ambig ./ max(max(ambig))))
view (-30,55);
axis tight
ylabel('\bf frequency')
xlabel('\bf delay')
zlabel('\bf ambiguity function')
figure
Nhalf = (size(ambig,1)-1)/2;
plot(delay,ambig(Nhalf+1,:)/(max(max(ambig))),'k')
xlabel('\bf delay')
ylabel('\bf normalized ambiguity cut for f=0')
grid
axis tight
figure
contour(delay,freq,(ambig ./ max(max(ambig))))
axis tight
ylabel('\bf frequency')
xlabel('\bf delay')
grid

end

