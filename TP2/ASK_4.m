clear all; close all;clc;

%Etudes de chaines de transmission     4_ASK     

Nb = 12000;
M = 4;
alpha = 0.5;
Fe = 12*10^3;
Rb = 6 * 10^3;
Rs = Rb / log2(M) ;
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
Symboles = symboles_decimaux.';
Symboles( symboles_decimaux == 0 ) = -3;
Symboles( symboles_decimaux == 1 ) = -1;
Symboles( symboles_decimaux == 3 ) = 1;
Symboles( symboles_decimaux == 2 ) = 3;


% Tracé des constellations en sortie du mapping 
figure();scatterplot(Symboles,1,0,'ro');

%suite des impultions dirac
Suite_dirac = kron(Symboles,[1 zeros(1,Ns-1)]);

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
z_echan = z(1:Ns:Ns*Nb/2);

%décisions
z_dec = z_echan;
z_dec(z_echan >= 2) = 3;
z_dec(z_echan <= -2) = -3;
z_dec(0 < z_echan & z_echan< 2) = 1;
z_dec(-2 < z_echan & z_echan <= 0) = -1;

%demapping
% on a besoin de faire le demapping car Teb = Tes / log2(M) (grace a gray)

%calcul du TES
%l'énergie en sortie du filtre de réception n'est pas maximale car la convolution porte sur moins de symboles.
%Du coup, quand les premiers symboles sont des +/- 3, ça génère une erreur Pour ne pas le prendre en compte 
%dans le calcul du TES je  les calculer à partir du 2e ou du 3e symbole

TES4 = sum(z_dec(2:end) ~= Symboles(2:end))/(length(Symboles) -1);
assert(TES4 == 0);

%%Implantation de la chaine avec bruit
TEB_bruit = zeros(1,7);
k = 0;

EbN0db = [0:6];
EbN0 = 10.^(EbN0db/10) ;

for k=1:length(EbN0) 
    
    %generation des bits
    bits = randi([0 1], 1,Nb);
    
    % groupement des bits par log2(M)
    blocs_bits=reshape(bits,log2(M),length(bits)/log2(M)).';

    % conversion des groupes de bits en une valeur décimale (de 0 à M-1)
    symboles_decimaux=bi2de(blocs_bits);

    %mapping
    Symboles = symboles_decimaux.';
    Symboles( symboles_decimaux == 0 ) = -3;
    Symboles( symboles_decimaux == 1 ) = -1;
    Symboles( symboles_decimaux == 3 ) = 1;
    Symboles( symboles_decimaux == 2 ) = 3;
    
    
    %suite des impultions dirac
    Suite_dirac = kron(Symboles,[1 zeros(1,Ns-1)]);
    
    %filtrage mise en forme
    xe = conv(Suite_dirac,h,'same');
        
    %calcul de la puissance du signal transmis
    Pr = mean(abs(xe).^2) ;
    
    %calcul du signal bruit 
    n = (sqrt((Pr*Ns)/(2*log2(M)*EbN0(k)))*randn(1,length(xe)));
        
    %ajout du bruit 
    xe = xe + n ;
    
    %filtrage de reception
    z =  conv(xe,h,'same');

    %echantillionage
    z_echan = z(1:Ns:Ns*Nb/2);
    
    %tracé des constellations
    figure(1);
    subplot(3,3,k+1);
    plot(real(z_echan),imag(z_echan),'k.');
    xlabel('I')
    ylabel('Q')
    title(["constellations pour Eb/No = ",k-1,"db"])

    %décisions
    z_dec = z_echan;
    z_dec(z_echan >= 2) = 3;
    z_dec(z_echan <= -2) = -3;
    z_dec(0 < z_echan & z_echan< 2) = 1;
    z_dec(-2 < z_echan & z_echan <= 0) = -1;
    
    %calcul du teb
    TEB_bruit(k) = sum(z_dec(2:end) ~= Symboles(2:end))/(length(Symboles) -1) / 2;
    
end

%Tracé du TEB calculé
figure;
semilogy(EbN0db,TEB_bruit);
xlabel('Eb/N0 db'),
ylabel('TEB'),
title('Tracé du TEB calculé');

%Comparaison du TEB theorique et TEB calculé
TEB_theo = (2 / log2(M)) * ((M-1) / M) * qfunc(sqrt(((6*log2(M))/(M^2-1))* 10.^(EbN0db/10)));
figure;
semilogy(EbN0db,TEB_theo);hold on;
semilogy(EbN0db,TEB_bruit);
xlabel('Eb/N0 db'),
ylabel('TEB'),
legend('TEB theorique','TEB calculé')
title('Comparaison du TEB theorique et TEB calculé ');
