clear all; close all;clc;
%################################################################################
%#               Etudes de chaines de transmission en bande de base             #
%################################################################################
%% Deuxième chaine : impact du choix du filtre de réception

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
x2=filter(h1,1,Suite_dirac);

%filtre de reception
hr2 = [ones(1,Ns/2) zeros(1,Ns/2)];

%filtrage reception
z2 = filter(hr2,1,x2);

%tracé du signal sortie du filtre reception
figure; plot([0:Te:(Nb*Ns-1)*Te],z2);
axis([0 10*Ts -Ns Ns]);  %pour une meilleur visualisation du signal
ylabel('z2');
title('le signal sortie du filtre reception (chaine 2)');

%diagrame de l'oeil
eye = reshape(z2,2*Ns,[]);
figure;plot(eye,'LineWidth',3);
ylabel ('z2(t)');
title("diagrame de l'oeil - (chaine 2)");

%echantillionage
z2_echan = z2(Ns:Ns:Ns*Nb);

%décisions
z2_dec = z2_echan > 0;

%calcul du TEB
TEB2 = sum(bits~=z2_dec)/length(bits); 
assert(TEB2 == 0);

%Implantation de la chaine avec bruit
TEB_bruit2 = zeros(1,7);
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
    z_temp = filter(hr2,1,xtemp+(sqrt((Pr*Ns)/(2*10^(i/10))))*randn(1,length(xtemp)));
    %echantiollanage
    z_echan_temp = z_temp(Ns:Ns:Ns*Nb);
    %decision
    z_dec_temp = z_echan_temp > 0;
    %calcul de TEB
    err = sum(bits_temp~=z_dec_temp) ;  
        if err >= 50 
        TEB_bruit2(i+1) = err/length(bits_temp);
        i = i+1;
        end
end

%tracé du TEB calculé 
figure;
semilogy(0:6,TEB_bruit2);
xlabel('Eb/N0 db'),
ylabel('TEB'),
title('tracé du TEB calculé échelle log (chaine 2)');

%Comparaison du TEB theorique et TEB calculé
TEB_theo2 = qfunc(sqrt(10.^([0:6]/10)));
figure;
semilogy(0:6,TEB_theo2);hold on;
semilogy(0:6,TEB_bruit2);
xlabel('Eb/N0 db'),
ylabel('TEB'),
legend('TEB theorique','TEB calculé')
title('Comparaison du TEB theorique et TEB calculé (chaine 2)');

%Comparaison du TEB calculé de la chaine 2 avec le TEB theorique de la chaine 1
TEB_theo1 = qfunc(sqrt(2*10.^([0:6]/10)));
figure;
semilogy(0:6,TEB_theo1);hold on;
semilogy(0:6,TEB_bruit2);
xlabel('Eb/N0 db'),
ylabel('TEB'),
legend('TEB calculé chaine 1','TEB calculé chaine 2')
title('Comparaison du TEB calculé de la chaine 2 avec le TEB calculé de la chaine 1');

