clear all; close all;clc;
%################################################################################
%#               Etudes de chaines de transmission en bande de base             #
%################################################################################
%% Troisième chaine à étudier : impact du choix du filtre de mise en forme

Rs = 3000;
alpha = .5;
Fe = 12000;
Te = 1/Fe;
Ts = 1/Rs;
Ns = Ts/Te;
Nb = 10000 ;                 
span = 8 ;

%generation des bits
bits = randi([0 1], 1,Nb);

%mapping
Symboles =2*bits-1;

%suite des impultions dirac
Suite_dirac = kron(Symboles,[1 zeros(1,Ns-1)]);

%filtre mise en forme
h3 = rcosdesign(alpha,span,Ns,'sqrt');

%filtrage mise en forme
x3 = conv(Suite_dirac,h3,'same');

%filtrage reception
z3 = conv(x3,h3,'same');

%tracé du signal sortie du filtre reception
figure; plot(z3);
axis([0 70*Ns -Ns Ns]);  %pour une meilleur visualisation du signal
xlabel('t en s');
ylabel('z3(t)');
title('le signal sortie du filtre reception (chaine 3)');

%diagrame de l'oeil
eye = reshape(z3,2*Ns,[]);
figure;plot(eye);
ylabel ('z3(t)');
title("diagrame de l'oeil - (chaine 3)");

%echantillionage
z3_echan = z3(1:Ns:Ns*Nb);

%décisions
z3_dec = z3_echan > 0;

%calcul du TEB
TEB3 = sum(bits~=z3_dec)/length(bits); 
assert(TEB3 == 0);

%%Implantation de la chaine avec bruit
TEB_bruit3 = zeros(1,7);
i = 0;
while i < 7 
    %generation des bits 
    bits_temp = randi([0 1], 1,Nb);
    %surechantillonage
    dirac_temp = kron(2*bits_temp-1,[1 zeros(1,Ns-1)]);
    %mise en forme
    xtemp=conv(dirac_temp,h3,'same');
    %ajout du bruit
    Pr = mean(abs(xtemp).^2);
    z_temp = conv(xtemp+(sqrt((Pr*Ns)/(2*10^(i/10))))*randn(1,length(xtemp)),h3,'same');
    %echantiollanage
    z_echan_temp = z_temp(1:Ns:Ns*Nb);
    %decision
    z_dec_temp = z_echan_temp > 0;
    %calcul de TEB
    err = sum(bits_temp~=z_dec_temp) ;
    if err >= 20
        TEB_bruit3(i+1) = err/length(bits_temp);
        i = i+1;
    end
    
end

%tracé du TEB calculé 
figure;
semilogy(0:6,TEB_bruit3);
xlabel('Eb/N0 db'),
ylabel('TEB'),
title('tracé du TEB calculé échelle log (chaine 3)');

%Comparaison du TEB theorique et TEB calculé
TEB_theo3 = qfunc(sqrt(2*10.^([0:6]/10)));
figure;
semilogy(0:6,TEB_theo3);hold on;
semilogy(0:6,TEB_bruit3);
xlabel('Eb/N0 db'),
ylabel('TEB'),
legend('TEB theorique','TEB calculé');
title('Comparaison du TEB theorique et TEB calculé (chaine 3)');

%Comparaison du TEB (chaine 3) et TEB (chaine 1)
TEB_theo1 = qfunc(sqrt(2*10.^([0:6]/10)));
figure;
semilogy(0:6,TEB_theo1);hold on;
semilogy(0:6,TEB_bruit3);
xlabel('Eb/N0 db'),
ylabel('TEB'),
legend('TEB théorique chaine 1','TEB calculé chaine 3')
title('Comparaison du TEB calculé de la chaine 3 avec le TEB théorique de la chaine 1');

%calcul de densité spectrale
DSP_x3 = (1/length(x3))*abs(fft(x3,2^nextpow2(length(x3)))).^2;

h1 = ones(1,Ns);

x=filter(h1,1,Suite_dirac);

DSP_x = (1/length(x))*abs(fft(x,2^nextpow2(length(x)))).^2;

%comparaison des DSP
figure;
plot(linspace(-Fe/2,Fe/2,2^nextpow2(length(x))),fftshift(DSP_x));hold on;
plot(linspace(-Fe/2,Fe/2,2^nextpow2(length(x3))),fftshift(DSP_x3));
ylabel ('DSP du signal');
xlabel('f en HZ');
legend('chaine1','chaine3');
title('la densité spectrale de puissance du signal transmis (chaine 1) vs (chaine 3)');

%introduction d'un canal de transmition BW = 1500
BW = 1500;
N = 20 ;

%filtre passe bas
hc=2*(BW/Fe)*sinc(2*BW*[-N/Fe:1/Fe:N/Fe]);

%filtrage canal  
c = conv(x3,hc,'same');

%filtrage reception
z3 = conv(c,h3,'same');

%diagrame de l'oeil
eye = reshape(z3,2*Ns,[]);
figure;plot(eye);
ylabel ('z3(t)');
title("diagrame de l'oeil - (chaine 3) avec filtrage canal BW = 1500");

%introduction d'un canal de transmition BW = 3000
BW = 3000;
N = 20 ;

%filtre passe bas
hc=2*(BW/Fe)*sinc(2*BW*[-N/Fe:1/Fe:N/Fe]);

%filtrage canal  
c = conv(x3,hc,'same');

%filtrage reception
z3 = conv(c,h3,'same');

%diagrame de l'oeil
eye = reshape(z3,2*Ns,[]);
figure;plot(eye);
ylabel ('z3(t)');
title("diagrame de l'oeil - (chaine 3) avec filtrage canal BW = 3000");