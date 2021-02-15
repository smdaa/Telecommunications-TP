clear all; close all;clc;

%Etudes de chaines de transmission     Q_PSK equivalente     

Nb = 10000;
M = 4;
alpha = 0.35;
Fe = 12*10^3;
Rb = 6 * 10^3;
Rs = Rb / log2(M) ;Ts = 1/Rs;
Te = 1/Fe;
Ns = Ts/Te;
span = 8 ;

%generation des bits
bits = randi([0 1], 1,Nb);

%mapping
symboles = 2*bits(1:2:end)-1 + 1i * (2*bits(2:2:end)-1) ;

%suite des impultions dirac
Suite_dirac = kron(symboles,[1 zeros(1,Ns-1)]);

%filtre mise enforme + filtre reception
h = rcosdesign(alpha,span,Ns,'sqrt');

%filtrage mise en forme
xe = conv(Suite_dirac,h,'same');

%tracé du signal en phase 
figure; plot(real(xe));
axis([0 70*Ns -1 1]);  %pour une meilleur visualisation du signal
ylabel('I(t)');
title('tracé du signal en phase ');

%tracé du signal en quadrature
figure; plot(imag(xe));
axis([0 70*Ns -1 1]);  %pour une meilleur visualisation du signal
ylabel('Q(t)');
title('tracé du signal en quadrature ');

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
z_echan = z(1:Ns:Ns*Nb/2);

%décisions + demapping
z_dec = zeros(1,Nb);
z_dec(1:2:end) = real(z_echan) > 0;
z_dec(2:2:end) = imag(z_echan) > 0;

%calcul du teb
TEB = sum(bits~=z_dec)/length(bits); 
assert(TEB == 0);

%%Implantation de la chaine avec bruit
TEB_bruit = zeros(1,7);
k = 0;

EbN0db = [0:6];
EbN0 = 10.^(EbN0db/10) ;

for k=1:length(EbN0) 
    
    %generation des bits
    bits = randi([0 1], 1,Nb);
    
    %mapping
    symboles = 2*bits(1:2:end)-1 + 1i * (2*bits(2:2:end)-1);
    
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
    
    %filtrage reception
    z =  conv(xe,h,'same');
    
    %echantillionnage
    z_echan = z(1:Ns:end);
    
    %tracé des constellations
    figure(1);
    subplot(3,3,k+1);
    plot(real(z_echan),imag(z_echan),'k.');
    xlabel('I')
    ylabel('Q')
    title(["constellations pour Eb/No = ",k-1,"db"])
    
    % on enregistre la constelation pour EB/N0 = 4 db
    if  k == 3
       util = z_echan; 
    end
    
    %decisions
    z_dec = zeros(1,Nb);
    z_dec(1:2:end) = real(z_echan) > 0;
    z_dec(2:2:end) = imag(z_echan) > 0;
    
    %calcul du TEB
    TEB_bruit(k) = sum(bits~=z_dec) / length(bits);
    
end

%Tracé des constellations en sortie du mapping et en sortie de l’échantillonneur
figure;
plot(real(util),imag(util), 'k.') ;hold on;
plot(real(symboles),imag(symboles), 'r*')
xlabel('I')
ylabel('Q')
legend('en sortie de l’échantillonneur','en sortie du mapping');
title("Tracé des constellations en sortie du mapping et en sortie de l’échantillonneur (Eb/N0 = 4)")

%Tracé du TEB calculé
figure;
semilogy(EbN0db,TEB_bruit);
xlabel('Eb/N0 db'),
ylabel('TEB'),
title('Tracé du TEB calculé');

%Comparaison du TEB theorique et TEB calculé
TEB_theo = qfunc(sqrt(2*10.^(EbN0db/10)));
figure;
semilogy(EbN0db,TEB_theo);hold on;
semilogy(EbN0db,TEB_bruit);
xlabel('Eb/N0 db'),
ylabel('TEB'),
legend('TEB theorique','TEB calculé')
title('Comparaison du TEB theorique et TEB calculé ');
