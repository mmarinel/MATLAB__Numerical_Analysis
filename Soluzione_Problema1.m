%     author: Matteo Marinelli
%     date: 03/03/2021
%     version: 1
%     note: Script per la risoluzione del Problema test 1; 
%           questo codice è composto di due parti: 
%               1. una pensata per la risoluzione del problema in generale,
%                  in cui, cioè, si risolve il problema sui dati forniti  
%                  nella prima parte del presente script, qualunque siano
%                  i valori di questi dati (i dati, quindi, possono anche
%                  essere diversi da quelli forniti nel testo)
%               2. e una in cui appunto si definiscono i particolari dati
%                  su cui si desidera far risolvere il problema al presente
%                  script durante la sua prossima esecuzione 
%                  (FASE DI SETTAGGIO).

%                 1. FASE DI SETTAGGIO 
                    
                    
% IMPORTANTE!!!:        
% Scrivere la funzione, per successive immissioni, in modo tale che agisca
% puntualmente sui vettori (i.e.: quando hai una funzione "complessa",
% ovvero una funzione scritta come combinazione di funzioni più semplici
% mediante operatori aritmetici, assicurati di usare le operazioni in modo
% puntuale qualora la loro natura non fosse già puntuale di per sé; 
% ad esempio: se hai * usa .*) 
f = @(x)sqrt(x);  
% qui sopra ho una funzione "semplice", quindi il discorso di
% prima non si applica
                    

% Nodi d'interpolazione                  
vx = (0:8).^2/64; 

% estremo sinistro dell'intervallo [a,b]
a = 0; 

% estremo destro dell'intervallo [a,b]
b = 1;            

% numero di punti in [a,b] utilizzati per i grafici
m = 1000; 

% Vettore (RIGA) di punti in [a,b] in cui valutare la 
% differenza tra polinomio di interpolazione e funzione da
% interpolare.
ts = (0:20)/20; 

% passo di graduazione sull'asse delle x
dex = 0.1; 

% passo di graduazione sull'asse delle y
dey = 2.5;  
 

%                 2. Parte "Generale"           


% Calcolo delle valutazioni di f nei nodi d'interpolazione vx
vy = f(vx);

% a pts verrà assegnato il vettore che contiene le valutazioni nei punti
% in ts del polinomio (p) d'interpolazione di f sui nodi vx
pts = valutazionePolinomioInterpolazione(vx, vy, ts);

differenze = pts - f(ts); % Vettore risultato delle differenze 

% Trasformo il risultato (differenze) in un vettore 'colonna'
differenze = transpose(differenze);

differenze %#ok<NOPTS>

%       PLOTTING
x = linspace(a, b, m); % Griglia di punti in [a,b]

% vettore delle valutazioni della funzione f ( sqrt(x) per il problema test
% 1) nei punti 'campione' dell'intervallo [0,1]
y = f(x);           

% vettore delle valutazioni del polinomio iterpolante p nei punti 'campione'
% dell'intervallo [0, 1]
pol = valutazionePolinomioInterpolazione(vx, vy, x); 

figure;
hold on;

% setting del plot della funzione y = f(x)
hp = plot(x, y, '-', 'LineWidth', 2);
% decido di dare a questo plot una colorazione verde scuro
set(hp, 'color', [0, 0.6, 0.2] ); 

% setting del plot del polinomio d'interpolazione
hp = plot(x, pol, '--', 'LineWidth', 3);
set(hp, 'color', 'Blue');

% setting del plot dei punti d'interpolazione
hp = plot(vx, vy, 'r.');
set(hp, 'MarkerSize', 25); 

%   ABBELLIMENTO AUTOMATICO DELLA FIGURA PRODOTTA
title('Problema 1', 'color', 'r');
box on;
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');
set(gca,'TickLabelInterpreter','latex');

% Settaggio del font
set(gca, 'fontsize', 18); 

% Settaggio della legenda
legend({'$f(x)$', 'polinomio d''interpolazione di $f(x)$', 'punti d''interpolazione'});

legend('Interpreter', 'latex');
legend('Location', 'northwest');

% Graduazione e stile degli assi

% OSS: Se voglio un passo delta (ad esempio = 0.1) negli assi graduati, ho
% che il passo rappresenta l'ampiezza di ciascun sottointervallo delimitato
% tra i punti campione che su quell'asse prendo per formare la 'scala'.
% Ora, l'ampiezza di ciasun sottointervallo (con i sottointervalli tutti di
% ampiezza uniforme) è proprio (b-a)/n, con n il numero di sottointervalli
% formati dall'asse graduato. Quindi passo = (b-a)/n --> n = (b-a)/passo, e
% #punti = n + 1 = (b-a)/passo  + 1.
ticks_x = linspace(a, b, ( ((b-a)/dex) + 1) ); 
xticks(ticks_x);
inizio  = f(a);
massimo = ceil(max(max(pol, y)));
ticks_y = linspace(inizio,  massimo, ( ((massimo-inizio)/dey) + 1) ); 
yticks(ticks_y);
set(gca, 'ticklabelinterpreter', 'latex');


%       Function Utilizzate nello Script

function pts = valutazionePolinomioInterpolazione(vx, vy, vt)


%     author: Matteo Marinelli
%     date: 28/02/2021
%     version: 2

%     input
%       vx: vettore dei nodi d'interpolazione
%       vy: vettore delle valutazioni di una certa funzione f sui nodi vx
%       vt: vettore dei punti in cui valutare il polinomio p interpolante f
%       sui nodi vx
%       
%     Output
%       pts: vettore delle valutazioni di p nei punti vt  


%inizializzazioni
n = length(vx);
m = length(vt);
pts = zeros(1, m); %inizializzazione del vettore di output


% Parte1 dell'algoritmo: calcolo delle differenze divise 
% (indipendente dai punti t)
dds = differenzEdivise(vx, vy);

% Parte2 dell'algoritmo (dipendente dai punti t)
for j = 1:m
    t = vt(j);
    hs = dds(n);
    for i = n-1:-1:1
        hs = dds(i) + (t - vx(i))*hs;
    end
    pts(j) = hs;
end

end

function dds = differenzEdivise(vx, vy)

%     author: Matteo Marinelli
%     date: 05/03/2021
%     version: 3

%     Input
%       vx: vettore (RIGA) dei nodi d'interpolazione
%       vy: vettore (RIGA) delle valutazioni di una certa funzione f 
%           sui nodi vx
%       
%       
%      Output  
%       dds: vettore (RIGA) dei coefficienti del polinomio che interpola,
%            nella sua forma di Newton, la funzione f sui nodi vx



%               Flusso Eccezionale

% Le cardinalità "dell'insieme" dei nodi di interpolazione e delle valutazioni
% della funzione f in quei nodi devono ovviamente coincidere!
if(length(vx) ~= length(vy))
    err('La cardinalità dell''insieme dei nodi di interpolazione e dell''insieme delle valutazioni della funzione f in quei nodi devono coincidere!')




end 


%           Flusso Standard

% quanti sono i nodi di interpolazione
n = length(vx);

% matrice n*n inizializzata a '0', corrispondente alla tabella delle 
% differenze divise
tabella = zeros(n);  


%           Calcolo della tabella delle differenze divise

% inizializzazione prima colonna
tabella(1:n, 1) = transpose(vy); 

% costruzione del resto della tabella
for i = 1:n
    for j = 2:i    
        % in (i,j) l'ultimo punto è xi e il penultimo è xj-1, quindi se levo 
        % l'ultimo elemento vado in riga j-1 e colonna (j-2)+1 dato che in 
        % colonna j ho f[x1,...,xj-1, xi] per ogni indice di riga i.
        tabella(i, j) = (tabella(i, j-1) - tabella(j-1, (j-2)+1))/(vx(i) - vx(j-1)); 
        
        % la quantità assegnata sopra è = f[x1, ..., xj-1, xi], con j < i.
    end 
end 

%           Calcolo del vettore dei coefficienti della forma di Newton
dds = transpose(diag(tabella)); 
   
end 



    

 










