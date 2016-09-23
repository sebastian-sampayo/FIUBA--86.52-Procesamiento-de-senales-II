function zeroplot(x,n_fig,titulo, color)
  figure(n_fig),
  w = linspace(0,2*pi);
  margen = 1.05*max(abs(x));
  axis([-margen margen -margen margen]);
  line(sin(w),cos(w),'LineStyle',':','Color',[0 0 0]);
  hold on;
  
  % scatter(real(x),imag(x), color);
  plot(real(x),imag(x), 'o', 'Color', color);
  grid on;
  axis square;
  xlabel('Re(z)','Fontsize',12);
  ylabel('Im(z)','Fontsize',12);
  title(titulo,'Fontsize',12);
  box
end
