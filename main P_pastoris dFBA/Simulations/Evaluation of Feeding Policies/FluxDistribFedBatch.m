clear all
model = readCbModel('PP_iFS618.xml');
load('metMovieFedBatch_decresing_mu_0.1_0.07.mat')

BiomassID = findRxnIDs(model,'BIOMASS');
GlucoseID = findRxnIDs(model,'EX_glc(e)');
PyrID = findRxnIDs(model,'EX_pyr(e)');
EtOHID = findRxnIDs(model,'EX_etoh(e)');
ArabID = findRxnIDs(model,'EX_abt_D(e)');
CitID = findRxnIDs(model,'EX_cit(e)');

% Rescatar mu
t = fluxDistrib(end,:);
mu = fluxDistrib(BiomassID,:);
glc = fluxDistrib(GlucoseID,:);
pyr = fluxDistrib(PyrID,:);
etoh = fluxDistrib(EtOHID,:);
arab = fluxDistrib(ArabID,:);
cit = fluxDistrib(CitID,:);

% Dibujar mu decreciente
t_mu = 35:1:100;
mu_dec_teor = zeros(1,length(t_mu));

for i=1:length(t_mu)
    % mu_dec_teor(i) = 0.10;
    mu_dec_teor(i) = 0.07+0.03*exp(-(t_mu(i)-35)*0.07);
end
figure(5)
subplot(2,3,1)
plot(t,mu,'b.',t_mu,mu_dec_teor,'r-')
xlabel('Time [h]')
ylabel('\mu h^-^1')

subplot(2,3,2)
plot(t,glc,'r.')
title('Glucose')
xlabel('Time [h]')
ylabel('[mmol/g_D_C_Wh]')

subplot(2,3,3)
plot(t,pyr,'r.')
title('Pyruvate')
xlabel('Time [h]')
ylabel('[mmol/g_D_C_Wh]')

subplot(2,3,4)
title('Ethanol')
plot(t,etoh,'r.')
xlabel('Time [h]')
ylabel('[mmol/g_D_C_Wh]')

subplot(2,3,5)
title('Arabitol')
plot(t,arab,'r.')
xlabel('Time [h]')
ylabel('[mmol/g_D_C_Wh]')

subplot(2,3,6)
title('Citrate')
plot(t,cit,'r.')
xlabel('Time [h]')
ylabel('[mmol/g_D_C_Wh]')