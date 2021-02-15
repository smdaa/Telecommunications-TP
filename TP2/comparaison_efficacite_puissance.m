clear all; close all;clc;

M_4ASK = 4;
M_16QAM = 16;
M_8PSK = 8;

EbN0db = [0:6];
EbN0 = 10.^(EbN0db/10) ;

TEB_theo_8PSK = ( 2 / log2(M_8PSK) ) * qfunc(sqrt(2 * log2(M_8PSK) * 10.^(EbN0db/10)) .* sin(pi/M_8PSK));
TEB_theo_16QAM = (4 / log2(M_16QAM)) * (1 - (1/sqrt(M_16QAM))) * qfunc(sqrt((3*log2(M_16QAM)/(M_16QAM-1))* 10.^(EbN0db/10)));
TEB_theo_4ASK = (2 / log2(M_4ASK)) * ((M_4ASK-1) / M_4ASK) * qfunc(sqrt(((6*log2(M_4ASK))/(M_4ASK^2-1))* 10.^(EbN0db/10)));
TEB_theo_QPSK = qfunc(sqrt(2*10.^(EbN0db/10)));

figure;
grid on ;
semilogy(EbN0db,TEB_theo_4ASK,'-s','LineWidth',2);hold on;
semilogy(EbN0db,TEB_theo_QPSK,'-o','LineWidth',1);hold on;
semilogy(EbN0db,TEB_theo_8PSK,'-x','LineWidth',1);hold on;
semilogy(EbN0db,TEB_theo_16QAM,'-*','LineWidth',1);
grid on;
xlabel('Eb/N0 db'),
ylabel('TEB'),
legend('4ASK','QPSK', '8PSK', '16QAM');
title('Comparaison des TEB des différentes chaines de transmission ');


