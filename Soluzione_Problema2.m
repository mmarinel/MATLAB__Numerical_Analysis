%     author: Matteo Marinelli
%     date: 05/03/2021
%     version: 1
%     note: Script per la risoluzione del Problema test 2; 
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




%               FASE DI SETTAGGIO

% l'estremo sinistro dell'intervallo di integrazione
a = 0; 

% l'estremo destro   dell'intervallo di integrazione
b = 1;          
                
% Nomi desiderati per le colonne della tabella del punto (c)                
colNames = {'I2', 'I4', 'I8', 'I16', 'I'};

% Nomi desiderati per le colonne della tabella da costruire nel punto (d)                
colNames2 = {'I2', 'I4', 'I8', 'I16', 'I', 'p0'};


% vettore delle precisioni fornite dal
% punto (b) del testo.
epss = linspace(10,10,10).^(-(1:10)); 

% vettore degli ordini delle formule
% dei trapezi citate nel punto (c)
ns   = linspace(2,2,4).^(1:4);        

% Punto in cui calcolare il polinomio d'interpolazione del 
% punto (d) del Problema.
pt_ptd = 0; 
 

%               Parte 'Generale'

fprintf("\nRisoluzione del problema 2\n\n");
%-------------------------------------------------------------------------

% Punto (a)
fprintf("Risoluzione punto (a):\n");
fprintf("Il punto (a) serve per la risoluzione del punto (b).\nLa risoluzione del punto (a) consiste nella function 'nepsilon' compresa nel presente script\n\n"); 





%-------------------------------------------------------------------------

% Punto (b)
fprintf("Risoluzione punto (b):\n");

% Inizializzazioni

% Lunghezza del vettore delle precisioni fornite 
% dal problema.
m  = length(epss); 

                           

% Costruzione Tabella

% Valore esatto dell'integrale fornito dal problema
I = exp(b) - exp(a);

% La colonna (VETTORE) della tabella richiesta con i valori esatti 
% dell'integrale 
I = linspace(I, I, m); 

% Vettore dei n(eps) per ogni eps in epss vettore di precisioni desiderate.
n_epss = nepsilon(a, b, epss);

% Vettore delle approssimazioni di I mediante le formule dei trapezi di
% ordine n_eps per ogni n_eps in n_epss
I_ns  = formulaTrapezi(@(x)exp(x), a, b, n_epss);

% Vettore degli errori |I-In|
err_I_In = abs(I - I_ns);

% Colonne della tabella
I        = transpose(I);
n_epss   = transpose(n_epss);
I_ns     = transpose(I_ns);
err_I_In = transpose(err_I_In);
epss     = transpose(epss);


table(n_epss, I_ns, I, err_I_In, epss)
%-------------------------------------------------------------------------

% Punto (c)
fprintf("Risoluzione punto (c):\n");

I    = exp(b) - exp(a);
I_ns = formulaTrapezi(@(x)exp(x), a, b, ns);

I_ns(length(I_ns)+1) = I; % Aggiungo un elemento: I

% Creazione e visualizzazione della tabella
array2table(I_ns, 'VariableNames', colNames)
%-------------------------------------------------------------------------

% Punto (d)
fprintf("Risoluzione punto (d):\n");

% Questi sono i nodi d'interpolazione corrispondenti ai passi di
% discretizzazione al quadrato delle formule dei trapezi di ordine n per
% approssimare l'integrale definito di f in [a,b], per ogni n in ns.
% (per ns si guardi la parte di SETTAGGIO)
num_nodes = length(ns);

% vettore dei nodi d'interpolazione
% citati nel punto d
vx   = (linspace( (b-a),(b-a),num_nodes )./(ns)).^(2); 

% Vettore delle approssimazioni mediante
% formula dei trapezi di ordine n > 1. (Ho
% tolto il valore esatto I).
vy   =  I_ns(1:length(I_ns)-1); 
                                
p0   =  valutazionePolinomioInterpolazione(vx, vy, pt_ptd);

I_ns(length(I_ns)+1) = p0;

array2table(I_ns, 'VariableNames', colNames2)
%------------------------------------------------------------------------



%               FUNCTION UTILIZZATE

%------------------------------------------------------------------------
function nepss = nepsilon(a, b, tolls)
   
%     author: Matteo Marinelli
%     date: 05/03/2021
%     version: 2
%     note: Codice per la risoluzione del punto (a) del secondo problema di
%           test. Il codice è pensato per un uso generale in cui passo al
%           programma estremi e tolleranze di valore qualunque.

%     Input
%       a:  l'estremo sinistro di un intervallo contenuto in R
%       b:  l'estremo destro di un intervallo   contenuto in R 
%   tolls:  vettore delle soglie di precisione desiderate per 
%           l'approssimazione, mediante formula dei trapezi di ordine n generico, 
%           del valore esatto dell'integrale definito della funzione 
%           exp(x) in [a,b]. 
%       
%     Output
%       neps: un vettore di costanti neps in R t.c. 
%             per ogni eps in tolls |I - In| <= eps per ogni n >=
%             neps. Dove con I intendiamo il valore esatto dell'integrale
%             definito di exp(x) in [a,b] e con In intendiamo
%             l'approssimazione di questo mediante la formula dei trapezi
%             di ordine n.
%

% nel caso fornissi gli estremi dell'intervallo nell'ordine sbagliato...
if( a > b )
    t = a;
    a = b;
    b = t;
end 

% OSS: |D(D(exp(x)))| = exp(x) <= exp(b) per ogni x in [a,b]


C = sqrt( ( ((b-a)^3)*exp(b) )/12 ); % Parte indipendente da 'tolls'

% N.B: L'operatore puntuale './' fa sì che venga 
% restituito un vettore in neps
nepss = ceil(C./sqrt(tolls)); 

end


function I_ns = formulaTrapezi(f, a, b, ns)

%     author: Matteo Marinelli
%     date: 05/03/2021
%     version: 3
        
%     Input
%       a: l'estremo sinistro di un intervallo [a,b]
%       b: l'estremo destro   di un intervallo [a,b]
%       ns: vettore di numeri interi positivi, relativi al numero di
%           sottointervalli di equal misura in cui suddividere l'intervallo
%           [a,b]; ovvero relativi all'ordine di una formula dei trapezi.
%       f: una certa funzione integrabile nell'intervallo di estremi a,b
%     
%     Output:
%       I_ns: vettore delle approssimazioni dell'integrale definito di f 
%             nell'intervallo di estremi a, b, mediante le formule dei 
%             trapezi di ordine n, per ogni n in ns.


hs = (b-a)./ns; % Passi di discretizzazione (VETTORE)

nh = length(hs);
I_ns = zeros(1, nh);

c = (f(a) + f(b))/2;

for i = 1:nh
    sums = sum(f( a + (1:ns(i)-1)*hs(i) ));
    
    % Valore approsimato dalla formula dei trapezi di ordine ns(i)
    % dell'integrale assegnato
    I_ns(i) = hs(i)*(c + sums); 
end                               
    

end

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



    

 


