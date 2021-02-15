clear all; close all;clc;
%################################################################################
%#               Etudes de chaines de transmission en bande de base             #
%################################################################################
%% Quatrième chaine à étudier : impact du choix du mapping

Nb = 10000 ;                   
Ns = 8;                        
Fe = 12000;
Te = 1/Fe;
Ts = Te*Ns;

%generation des bits
bits = randi([0 1], 1,Nb);

%filtre mise en forme
h4 = ones(1,Ns);

%mapping
Symboles = (2 * bi2de(reshape(bits, 2, length(bits)/2).','left-msb') - 3).';

%suite des impultions dirac
Suite_dirac = kron(Symboles,[1 zeros(1,Ns-1)]);

%mise en forme
x4 = filter(h4,1,Suite_dirac);
figure; plot([0:Te:(Nb*Ns-1)*Te/2],x4);
axis([0 20*Ts -4 4]);  %pour une meilleur visualisation du signal
xlabel('t en s');
ylabel('x(t)');
title('le signal a transmettre apres le filtrage mise en forme (chaine 4)');

%calcul de densité spectrale
DSP_x4 = (1/length(x4))*abs(fft(x4,2^nextpow2(length(x4)))).^2;

%affichage de la DSP
figure; 
plot(linspace(-Fe/2,Fe/2,2^nextpow2(length(x4))),fftshift(DSP_x4));
ylabel ('DSP du signal');
xlabel('f en HZ');
title('la densité spectrale de puissance du signal transmis (chaine 4)');

%comparaison des DSP
x=filter(h4,1,kron(2*bits-1,[1 zeros(1,Ns-1)]));
DSP_x = (1/length(x))*abs(fft(x,2^nextpow2(length(x)))).^2;
figure;
plot(fftshift(DSP_x4));hold on;
plot(fftshift(DSP_x));
ylabel ('DSP du signal');
legend('chaine4','chaine1');
title('la densité spectrale de puissance du signal transmis (chaine 1) vs (chaine 4)');

%filtrage reception
z4 = filter(h4,1,x4);

%diagrame de l'oeil
eye = reshape(z4,2*Ns,[]);
figure;plot(eye,'LineWidth',3);
ylabel ('z4(t)');
title("diagrame de l'oeil - (chaine 4)");

%echantillionage
z4_echan = z4(Ns:Ns:Ns*Nb/2);

%décisions
z4_dec = z4_echan;
z4_dec(z4_echan >= 2*Ns) = 3;
z4_dec(z4_echan <= -2*Ns) = -3;
z4_dec(0 < z4_echan & z4_echan< 2*Ns) = 1;
z4_dec(-2*Ns < z4_echan & z4_echan<= 0) = -1;


%demapping
BitsDecides = reshape(de2bi((z4_dec + 3)/2,'left-msb').',1,length(bits));

%calcul du TEB
TEB4 = sum(bits~=BitsDecides)/length(bits); 
TES4 = sum(z4_dec ~= Symboles)/length(Symboles);
assert(TEB4 == 0);
assert(TES4 == 0);

    %Implantation de la chaine avec bruit - calcul de TES / TEB
    TEB_bruit4 = zeros(1,7);
    TES_bruit4 = zeros(1,7);
    i = 0;
    while i < 7
        %generation des bits 
        bits_temp = randi([0 1], 1,Nb);

        %mapping
        Symboles = (2 * bi2de(reshape(bits_temp, 2, length(bits_temp)/2).','left-msb' ) - 3).';

        %suite des impultions dirac
        Suite_dirac = kron(Symboles,[1 zeros(1,Ns-1)]);

        %mise en forme
        xtemp = filter(h4,1,Suite_dirac);

        %ajout du bruit
        Pr = mean(abs(xtemp).^2);

        %filtrage reception
        ztemp = filter(h4,1,xtemp+(sqrt((Pr*Ns)/(4*10^(i/10))))*randn(1,length(xtemp)));

        %echantillionage
        z4_echan_temp = ztemp(Ns:Ns:end);

        %décisions
        z4_dec = z4_echan_temp;
        z4_dec(z4_echan_temp >= 2*Ns) = 3;
        z4_dec(z4_echan_temp <= -2*Ns) = -3;
        z4_dec(0 < z4_echan_temp & z4_echan_temp< 2*Ns) = 1;
        z4_dec(-2*Ns < z4_echan_temp & z4_echan_temp<= 0) = -1;

        %demapping
        BitsDecides = reshape(de2bi((z4_dec + 3)/2,'left-msb').',1,length(bits_temp));
        errb = sum(bits_temp~=BitsDecides);
        errs = sum(z4_dec~=Symboles);
            if errs > 100
            TES_bruit4(i+1) = errs/length(Symboles);
            TEB_bruit4(i+1) = errb/length(bits_temp);
            i = i+1;
            end
    end

%tracé du TES calculé 
figure;
semilogy(0:6,TES_bruit4);
xlabel('Eb/N0 db');
ylabel('TES');
title('tracé du TES calculé  (chaine 4)');

%tracé du TEB calculé 
figure;
semilogy(0:6,TEB_bruit4);
xlabel('Eb/N0 db');
ylabel('TEB');
title('tracé du TEB calculé  (chaine 4)');

%Comparaison du TES chaine theorique et TES calculé
TES_theo4 = (3/2)*qfunc(sqrt((4/5)*10.^([0:6]/10)));
figure;
semilogy(0:6,TES_theo4);hold on;
semilogy(0:6,TES_bruit4);
xlabel('Eb/N0 db');
ylabel('TES');
legend('TES theorique','TES calculé')
title('Comparaison du TES theorique et TES calculé (chaine 4)');

%Comparaison du TEB chaine theorique et TEB calculé
TEB_theo4 = TES_theo4/2;  %cas de gray dans la partie theorique
figure;
semilogy(0:6,TEB_theo4);hold on;
semilogy(0:6,TEB_bruit4);
xlabel('Eb/N0 db');
ylabel('TES');
legend('TEB theorique','TEB calculé')
title('Comparaison du TEB theorique et TEB calculé (chaine 4)');

