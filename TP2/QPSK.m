clear all; close all;clc;

%Etudes de chaines de transmission     Q_PSK   

Nb = 10000;
alpha = 0.35;
fp = 2*10^3;
Fe = 10*10^3;
Rs = 10^3;
Ts = 1/Rs;
Te = 1/Fe;
Ns = Ts/Te;
span = 8 ;
M = 4;
T = [0:Te:(Ns*Nb/2 -1)*Te];

%generation des bits
bits = randi([0 1], 1,Nb);

%mapping
symboles = 2*bits(1:2:end)-1 + 1i * (2*bits(2:2:end)-1);

% Tracé des constellations en sortie du mapping 
scatterplot(symboles,1,0,'ro');

%suite des impultions dirac
Suite_dirac_I = kron(real(symboles),[1 zeros(1,Ns-1)]);
Suite_dirac_Q = kron(imag(symboles),[1 zeros(1,Ns-1)]);

%filtre mise enforme + filtre reception
h = rcosdesign(alpha,span,Ns,'sqrt');

%filtrage mise en forme
I = conv(Suite_dirac_I,h,'same');
Q = conv(Suite_dirac_Q,h,'same');

%signal transmis sur fréquence porteuse
xe = I+1i*Q;

%transposition de fréquence
xe = xe.*exp(2*1i*pi*fp*T);

%tracé du signal en phase 
figure; plot(I);
axis([0 70*Ns -1 1]);  %pour une meilleur visualisation du signal
ylabel('I(t)');
title('tracé du signal en phase ');

%tracé du signal en quadrature
figure; plot(Q);
axis([0 70*Ns -1 1]);  %pour une meilleur visualisation du signal
ylabel('Q(t)');
title('tracé du signal en quadrature ');

%tracé du signal transmis sur fréquence porteuse
figure; plot(real(xe));
axis([0 20*Ns -1 1]);  %pour une meilleur visualisation du signal
ylabel('xe');
title('tracé du signal transmis sur fréquence porteuse');

% DSP du signal à transmettre
DSP = (1/length(real(xe)))*abs(fft(real(xe),2^nextpow2(length(real(xe))))).^2;

%affichage de la DSP
figure; 
plot(linspace(-Fe/2,Fe/2,2^nextpow2(length(real(xe)))),fftshift(DSP));
ylabel ('DSP du signal');
xlabel('f en HZ');
title('la densité spectrale de puissance du signal transmis ');

%retour en bande de base
z_I = real(xe).*(cos(2*pi*fp*T));
z_Q = real(xe).*(sin(2*pi*fp*T));

%filtrage passe bas 
hb = 2*(fp/Fe)*sinc(2*(fp/Fe)*[-20:20]); %filtre d'ordre 20 

z_I = conv(z_I,hb,'same');
z_Q = conv(z_Q,hb,'same');
z = z_I -1i * z_Q;

%filtrage de reception
z =  conv(z,h,'same');

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
    
    %transposition de fréquence
    xe = xe.*exp(2*1i*pi*fp*[0:Te:(Nb*0.5*Ns-1)*Te]);
    
    xe = real(xe);
    
    %calcul de la puissance du signal transmis 
    Pr = mean(abs(xe).^2);
    
    %ajout du bruit
    bruit = sqrt((Pr*Ns)/(2*log2(M)*EbN0(k)))*randn(1,length(xe));
    xe = xe + bruit;
    
    %retour en bande de base
    z_I = xe.*(cos(2*pi*fp*T));  
    z_Q = xe.*(sin(2*pi*fp*T));
    
    %filtrage passe bas
    z_I = conv(z_I,hb,'same');
    z_Q = conv(z_Q,hb,'same');

	z = z_I -1i * z_Q;
    
    %filtrage reception
    z =  conv(z,h,'same');
    
    %echantillionnage
    z_echan = z(1:Ns:end);
    
    %decisions
    z_dec = zeros(1,Nb);
    z_dec(1:2:end) = real(z_echan) > 0;
    z_dec(2:2:end) = imag(z_echan) > 0;
    
    %calcul TEB
    TEB_bruit(k) = sum(bits~=z_dec) / length(bits);
    
end

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
