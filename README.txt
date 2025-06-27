L'assignment è stato svolto articolando il procedimento in 4 declinazioni:
1) Studio del processo con drift di tipo k
2) Studio del processo con drift di tipo alpha
3) Studio del processo con drift di tipo alpha regolarizzato
4) Proof of concept, studio della possibilità di parziale parallelizzazione del calcolo

Le diverse cartelle contengono il codice in C (con uso della libreria math) nella specifica declinazione
e il file di input adatto. L'ottimizzazione della memoria non è massima, non si consiglia di generare
serie temporali più lunghe di 10^11 punti (max_temp*intervallini per le serie temporali,
 max_temp_ens*intevallini*numero serie per ensemble) a meno di disporre di una ram superiore a 32gb.
(i valori adatti potrebbero richiedere un minuto di fine tuning)

NB: per la parallelizzazione il codice impiega ~1gb di ram per processore per un tempo fisico di 
10^6 con 100dt, per 1024 serie dovrebbe completare il processo in 14min su 8 processori.
Compilato con openMPI 5.0.6.

Gabriele Sano
