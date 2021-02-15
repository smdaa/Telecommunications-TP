clear all; close all;clc;
%################################################################################
%#               Etudes de chaines de transmission en bande de base             #
%################################################################################
%% Premiere chaine: chaine de réference

Nb = 10000 ;                   
Ns = 8;                        
Fe = 12000;
Te = 1/Fe;
Ts = Te*Ns;

%generation des bits
bits = randi([0 1], 1,Nb);

%mapping
Symboles =2*bits-1;

%suite des impultions dirac
Suite_dirac = kron(Symboles,[1 zeros(1,Ns-1)]);

%filtre de mise en forme 
h1 = ones(1,Ns);

%filtrage mise en forme 
x=filter(h1,1,Suite_dirac);
figure; plot([0:Te:(Nb*Ns-1)*Te],x);
axis([0 7*Ts -1.5 1.5]);  %pour une meilleur visualisation du signal
xlabel('t en s');
ylabel('x(t)');
title('le signal a transmettre apres le filtrage mise en forme (chaine 1)');

%calcul de densité spectrale
DSP_x = (1/length(x))*abs(fft(x,2^nextpow2(length(x)))).^2;

%affichage de la DSP
figure; 
plot(linspace(-Fe/2,Fe/2,2^nextpow2(length(x))),fftshift(DSP_x));
ylabel ('DSP du signal');
xlabel('f en HZ');
title('la densité spectrale de puissance du signal transmis (chaine 1)');

%filtrage reception 
z1 = filter(h1,1,x);

%tracé du signal sortie du filtre reception
figure;
plot([0:Te:(Nb*Ns-1)*Te],z1);
axis([0 7*Ts -2*Ns 2*Ns]);  %pour une meilleur visualisation du signal
xlabel('t en s');
ylabel('z(t)');
title('le signal sortie du filtre reception (chaine 1)');

%diagrame de l'oeil
eye = reshape(z1,2*Ns,[]);
figure;plot(eye,'LineWidth',3);
ylabel ('z(t)');
title("diagrame de l'oeil - (chaine 1)");

%echantillionage
z1_echan = z1(Ns:Ns:Ns*Nb);

%décisions + demapping
z1_dec = z1_echan > 0;

%calcul du TEB
TEB1 = sum(bits~=z1_dec)/length(bits); 
assert(TEB1 == 0);  

%%Implantation de la chaine avec bruit
TEB_bruit1 = zeros(1,7);
i = 0;
while i < 7 
    %generation des bits 
    bits_temp = randi([0 1], 1,Nb);
    %surechantillonage
    dirac_temp = kron(2*bits_temp-1,[1 zeros(1,Ns-1)]);
    %mise en forme
    xtemp=filter(h1,1,dirac_temp);
    %ajout du bruit
    Pr = mean(abs(xtemp).^2);
    z_temp = filter(h1,1,xtemp+(sqrt((Pr*Ns)/(2*10^(i/10))))*randn(1,length(xtemp)));
    %echantiollanage
    z_echan_temp = z_temp(Ns:Ns:Ns*Nb);
    %decision
    z_dec_temp = z_echan_temp > 0;
    %calcul de TEB
    err = sum(bits_temp~=z_dec_temp) ;
    if err >= 20
        TEB_bruit1(i+1) = err/length(bits_temp);
        i = i+1;
    end
end

%tracé du TEB calculé 
figure;
semilogy(0:6,TEB_bruit1);
xlabel('Eb/N0 db'),
ylabel('TEB'),
title('tracé du TEB calculé échelle log (chaine 1)');

%Comparaison du TEB theorique et TEB calculé
TEB_theo1 = qfunc(sqrt(2*10.^([0:6]/10)));
figure;
semilogy(0:6,TEB_theo1);hold on;
semilogy(0:6,TEB_bruit1);
xlabel('Eb/N0 db'),
ylabel('TEB'),
legend('TEB theorique','TEB calculé')
title('Comparaison du TEB theorique et TEB calculé (chaine 1)');
