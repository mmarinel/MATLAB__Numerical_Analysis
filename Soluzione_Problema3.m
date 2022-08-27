%     author: Matteo Marinelli
%     date: 05/03/2021
%     version: 1
%     note: Script per la risoluzione del Problema test 3; 
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


%           FASE DI SETTAGGIO


% Sistema Lineare Ax = b da risolvere
A = [5,1,2; -1,7,1; 0,1,-3];
b = [13; 16; -7];

% Dati per il metodo di Jacobi 

% Vettore d'innesco
x0 = [0; 0; 0];

% Le prime n_itera iterazioni che si vogliono far fare al
% metodo nel punto (a).
n_itera = 10;   
 
% Vettore delle precisioni citate nel punto (c) della 
% traccia del problema test 3.
epss = 10.^(-(1:10)); 

% la norma che si vuole utilizzare nel punto (c)
norma = Inf; 

% I nomi desiderati per le colonne della tabella da costruire per la
% risoluzione del punto (c) del problema test 3.
ColNamesc = {'Keps_s', 'xe_s', 'xx_s', 'norm_Inf_err_s', 'eps_s__AKA_precisioni'};

%----------------------------------------------------------------


%           PARTE 'GENERALE'

% Punti (a) & (b)
fprintf("\nRisoluzione dei punti (a) & (b):\n")

x = A\b; % soluzione esatta del sistema

[~, n] = size(A);

% inizializzazione della matrice che conterrà le prime n_itera iterazioni
% del metodo di Jacobi
xs = zeros(n, n_itera); 
                  

% Calcolo delle prime n = n_itera iterazioni del metodo di Jacobi
% ogni iterazione andrà ad occupare una colonna nella matrice xs
for i = 1:n_itera
    [xs(:, i), ~, ~]  = jacobi(A, b, 0, x0, 2, i); 
    % Se metto toll = 0, allora avrò che o esco prima, ma in questo caso
    % soluzione esatta del sistema e quella calcolata al passo j < i
    % coincideranno, dunque anche se vado avanti fino all'iterazione i la
    % soluzione calcolata resterà la stessa ( perchè se il sistema converge 
    % allora è per forza consistente, e dalla consistenza ricavo l'osservazione
    % che le soluzioni calcolate resteranno invariate nelle successive
    % iterazioni del metodo), oppure esco esattamente all'iterazione i-esima
    % perchè i è il numero massimo di iterazioni consentite al metodo di
    % Jacobi e non convergo esattamente (cioè con tolleranza nulla) 
    % prima di i iterazioni 
    % (la soluzione potrebbe convergere esattamente anche all'iterazione i, 
    %  in ogni caso oltre i iterazioni non posso andare).
    % Si nota, intoltre, che usare la norma 2, piuttosto che la norma Inf o
    % la norma 1, è indifferente quando invoco la function che implementa
    % il metodo di Jacobi su una tolleranza di valore nullo. 
    % Questo perché, dalla proprietà di positività delle norme, si ha che 
    % la norma di un vettore è sempre >= 0, e inoltre è 0 sse il vettore 
    % stesso è 0, indipendentemente dalla particolare norma utilizzata.
end 

tabella = [x0, xs, x];

tabella %#ok<NOPTS> 
%-------------------------------------------------------------------------

% Punto (c)
fprintf("Risoluzione punto (c):\n")

% inizializzazioni 0
[~, n] = size(A);    % Quante componenti avranno i vettori soluzione del sistema
rows = length(epss); % Quante righe ci saranno nella tabella

% inizializzazioni 1
Kepss = zeros(rows, 1); % Vettore colonna dei Keps per ogni eps in epss
xes   = zeros(rows, n); % matrice in cui ogni riga corrisponde ad una soluzione 
                        % approssimata xeps calcolata dal metodo di Jacobi, 
                        % con eps in epss la precisione relativa a xeps
xx    = zeros(rows, n); % ogni riga corrisponde alla soluzione esatta

norm_err = zeros(rows, 1); % Ogni riga corrisponde alla norma 'norma' 
                               % (si guardi la fase di SETTAGGIO)
                               % della differenza tra la soluzione esatta x 
                               % e la soluzione approssimata xeps (relativa
                               % alla precisione eps in epss).

x = A\b; % Soluzione esatta del sistema (La ripeto per leggibilita')

%       COSTRUZIONE TABELLA


% Se voglio raggiungere una precisione del tipo
% norm( (x - x^(k)), 'norma') <= eps, allora ho che la soglia di
% arresto per il criterio d'arresto del residuo, rispetto al 
% condizionamento di A in norma 'norma', è
% eps/( norm(x, 'norma')*condizionamento(A, 'norma') ). Infatti, 
% norm(r^(k), 'norma')/norm(b, 'norma') <= toll, con toll la soglia di  
% arresto del criterio del residuo,  ==> (IMPLICA)
% norm(x - x^(k), 'norma')/ norm(x, 'norma') <= condizionamento(A, 'nor
% ma')*toll.
% Se pongo toll = eps/( norm(x,'norma')*condizionamento(A, 'norma') ),  
% allora ho che norm(x-x^(k), 'norma') <= norm(x,'norma')*condizioname 
% nto(A, 'norma')*toll = eps.
tolls = epss/( norm(x, norma)*condizionamento(A, norma) );

for i = 1:rows
% So che mettere 'Inf' ( quarto parametro di jacobi) è un rischio, ma noi
% stiamo lavorando nell'ipotesi che il metodo converga con precisione
% 'tolls(i)' alla soluzione desiderata in un numero "ragionevole" di passi;
% altrimenti non potrei soddisfare la richiesta del punto (c), problema 3.
    [xt, k, ~] = jacobi(A, b, tolls(i), x0, norma, Inf); 
    
    Kepss(i, 1) = k;
    
    % Perchè xt è un vettore 'colonna' e a me serve un vettore 'riga'
    xes(i, :) = transpose(xt); 
    
    % Stessa ragione di su
    xx(i, :)   = transpose(x); 
    
    norm_err(i, 1) = norm(x-xt, norma);
end 

precisioni = transpose(epss);

tabella = table(Kepss, xes, xx, norm_err, precisioni, 'VariableNames', ColNamesc);

tabella %#ok<NOPTS>
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------

%...       FUNCTION UTILIZZATE
%-------------------------------------------------------------------------

function [xk, k, rk_norm] = jacobi(A, b, toll, x0, norma, Nmax)

%     author: Matteo Marinelli
%     date: 28/02/2021
%     version: 2

%     Input
%       A: una matrice quadrata
%       b: un vettore di termini noti relativo al sistema lineare da
%          risolvere Ax = b
%       toll: una soglia di precisione relativa alla condizione d'arresto
%             della function rispetto al criterio d'arresto del residuo con
%             tolleranza, appunto, 'toll'
%       x0:   un vettore di innesco
%       Nmax: un numero massimo di iterazioni consentite al programma
%      norma: la norma rispetto alla quale si considera il criterio
%             d'arresto del residuo.
%       
%     Output
%       xk :  primo vettore (della successione di vettori generata) che
%             soddisfa il criterio di arresto del residuo relativamente
%             alla tolleranza 'toll' e norma 'norma', oppure, 
%             l'ultimo vettore calcolato dal programma: (x^Nmax), qualora 
%             nessun vettore generato dovesse soddisfare il criterio di 
%             arresto sopra citato
%       k :   il numero di iterazioni eseguite dal programma
%  rk_norm:   la norma 'norma' dell'ultimo residuo calcolato

xk = x0;
rk = b - A*xk;
rk_norm = norm(rk, norma);  
b_norm = norm(b, norma);         
norm_rb = rk_norm/b_norm;
D = diag(diag(A)); % Precondizionatore

k = 0;
while( norm_rb > toll && k < Nmax)
    zk = D\rk;
    xk = xk + zk;
    rk = b - A*xk;
    rk_norm = norm(rk, norma);
    norm_rb = rk_norm/b_norm;
    k = k + 1;
end

end


function cnd = condizionamento(A, norma)

%    Input
%       A: una matrice quadrata
%       
%    Output
%       cnd: il condizionamento di A in norma 'norma'

cnd = norm(A, norma)*norm(inv(A), norma);

end
                     
                     
                     
                     