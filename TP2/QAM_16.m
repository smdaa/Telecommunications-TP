clear all; close all;clc;

%Etudes de chaines de transmission     16_QAM     

Nb = 16000;
M = 16;
alpha = 0.5;
Fe = 12*10^3;
Rb = 6 * 10^3;
Rs = Rb / log2(M);
Ts = 1/Rs;
Te = 1/Fe;
Ns = Ts/Te;
span = 8 ;


%generation des bits
bits = randi([0 1], 1,Nb);

% groupement des bits par log2(M)
blocs_bits=reshape(bits,log2(M),length(bits)/log2(M)).';

% conversion des groupes de bits en une valeur décimale (de 0 à M-1)
symboles_decimaux=bi2de(blocs_bits);

%mapping
symboles = qammod(symboles_decimaux,M,'gray').';

% Tracé des constellations en sortie du mapping 
figure();scatterplot(symboles,1,0,'ro');

%suite des impultions dirac
Suite_dirac = kron(symboles,[1 zeros(1,Ns-1)]);

%filtre mise enforme + filtre reception
h = rcosdesign(alpha,span,Ns,'sqrt');

%filtrage mise en forme
xe = conv(Suite_dirac,h,'same');

%DSP du l'envloppe complexe
DSP = (1/length(xe))*abs(fft(xe,2^nextpow2(length(xe)))).^2;
figure; 
plot(linspace(-Fe/2,Fe/2,2^nextpow2(length(real(xe)))),fftshift(DSP));
ylabel ('DSP du signal');
xlabel('f en HZ');
title('DSP du signal modulé sur porteuse');

%filtrage de reception
z =  conv(xe,h,'same');

%echantillionage
z_echan = z(1:Ns:end);

%décisions + demapping
z_dec = qamdemod(z_echan.',M,'gray');
demodbi=de2bi(z_dec.');
bits_decides=reshape(demodbi.',size(demodbi,1)*size(demodbi,2),1).';

TEB = sum(bits(3:end)~=bits_decides(3:end))/length(bits);

%%Implantation de la chaine avec bruit
TEB_bruit = zeros(1,7);
k = 0;

EbN0db = [0:6];
EbN0 = 10.^(EbN0db/10) ;

for k=1:length(EbN0) 
    
    % groupement des bits par log2(M)
    blocs_bits=reshape(bits,log2(M),length(bits)/log2(M)).';

    % conversion des groupes de bits en une valeur décimale (de 0 à M-1)
    symboles_decimaux=bi2de(blocs_bits);

    %mapping
    symboles = qammod(symboles_decimaux,M,'gray').';
    
    %suite des impultions dirac
    Suite_dirac = kron(symboles,[1 zeros(1,Ns-1)]);
    
    %filtrage mise en forme
    xe = conv(Suite_dirac,h,'same');
        
    %calcul de la puissance du signal transmis
    Pr = mean(abs(xe).^2) ;
    
    %calcul du signal bruit sur la voie I 
    n_I = (sqrt((Pr*Ns)/(2*log2(M)*EbN0(k)))*randn(1,length(xe)));
    
    %calcul du signal bruit sur la voie Q 
    n_Q = (sqrt((Pr*Ns)/(2*log2(M)*EbN0(k)))*randn(1,length(xe)));
    
    %ajout du bruit 
    xe = xe + (n_I + (1i * n_Q));
    
    
    %filtrage de reception
    z =  conv(xe,h,'same');

    %echantillionage
    z_echan = z(1:Ns:end);
    
    %tracé des constellations
    figure(1);
    subplot(3,3,k+1);
    plot(real(z_echan),imag(z_echan),'k.');
    xlabel('I')
    ylabel('Q')
    title(["constellations pour Eb/No = ",k-1,"db"])

    %décisions + demapping
    z_dec = qamdemod(z_echan.',M,'gray');
    demodbi=de2bi(z_dec.');
    bits_decides=reshape(demodbi.',size(demodbi,1)*size(demodbi,2),1).';

    
    %calcul du teb
    TEB_bruit(k) = sum(bits(3:end) ~= bits_decides(3:end))/(length(bits) -1) ; 
    
end

%Tracé du TEB calculé
figure;
semilogy(EbN0db,TEB_bruit);
xlabel('Eb/N0 db'),
ylabel('TEB'),
title('Tracé du TEB calculé');

%Comparaison du TEB theorique et TEB calculé
TEB_theo = (4 / log2(M)) * (1 - (1/sqrt(M))) * qfunc(sqrt((3*log2(M)/(M-1))* 10.^(EbN0db/10)));
figure;
semilogy(EbN0db,TEB_theo);hold on;
semilogy(EbN0db,TEB_bruit);
xlabel('Eb/N0 db'),
ylabel('TEB'),
legend('TEB theorique','TEB calculé')
title('Comparaison du TEB theorique et TEB calculé ');
