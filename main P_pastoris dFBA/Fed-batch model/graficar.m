function Graficar(t,x,odeEq,tFinal,Exp,labels0)

if ~exist(strjoin({'Resultados ',Exp},''), 'dir')
    mkdir(strjoin({'Resultados ',Exp},''));
end

warning(strjoin({'Los resultados se guardar�n en ',pwd,'\Resultados ',Exp,'\'},''));

for i = 1:length(odeEq)
    figure('Visible','off')
    fprintf('\nGraficando ');
    fprintf(strjoin(labels0(odeEq(i)+1)));
    
    for j = 1:length(tFinal)
        subplot(2,2,j)
        plot(t(1:ceil(length(t)/tFinal(j)),1),x(1:ceil(length(t)/tFinal(j)),odeEq(i)),'-');
        title(strjoin([Exp,labels0(odeEq(i)+1),'t =',num2str(ceil(t(ceil(length(t)/tFinal(j)),1))),'s']));
        axis tight
        orient landscape
        print('-dpdf','-r300',[pwd strjoin(['\Resultados ',Exp,'\','TerpenoidSim - ',num2str(i,'%02.0f'),' - ',labels0(odeEq(i)+1),'.pdf'],'')])
    end
end
