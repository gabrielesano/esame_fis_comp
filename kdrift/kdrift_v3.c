#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <limits.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// costanti
const double K = 1.0;
const double D = 1.0;

// generazione vettore
struct Vector {
    double *data;
    size_t size;
};

struct Vector* newVector(size_t sz) {
    struct Vector *retVal = malloc(sizeof(struct Vector));
    if (retVal == NULL) return NULL;
    retVal->data = malloc(sz * sizeof(double));
    if (retVal->data == NULL) { free(retVal); return NULL; }
    retVal->size = sz;
    return retVal;
}

void delVector(struct Vector *vector) {
    if (vector != NULL) {
        free(vector->data);
        free(vector);
    }
}

// prototipi funzioni
double box_muller();
void calcolo_stampa_momenti(struct Vector* vec, const char* series_type, FILE *outfile);
double max(struct Vector * array);
double min(struct Vector * array);
void istogramma(struct Vector * array, int binNum, const char *name);
void run_long_timeseries_analysis(long long int n_steps, double dt, double x0, int seed, int n_bins, int max_lag_seconds);
void calcolo_autocorr_ts(struct Vector* vec, int max_lag_seconds, const char* filename);
void calcolo_autocorr_ens_average(struct Vector** all_series, int n_series, int maxTime_ensemble, int intervallini, const char* filename);


/////////////////////////////////////////////////////////////////main////////////////////////////////////////////////////////////////
int main() {
    // dichiarazioni e lettura da file
    int maxTime_long, maxTime_ensemble, intervallini, N_series;
    int seed_ensemble, n_bins, seed_long, max_lag_seconds_from_file, ensemble_autocorr_toggle;
    double x0;


    FILE *infile = fopen("input_secondo.dat", "r");
    if (infile == NULL) {
        perror("Err: non esiste 'input_secondo.dat'");
        return 1;
    }

    // Ordine parametri input: maxTime_long maxTime_ensemble intervallini N_series n_bins x0 seed_ensemble seed_long max_lag_seconds ensemble_autocorr_toggle
    int numeri_input = fscanf(infile, "%d %d %d %d %d %lf %d %d %d %d",
                            &maxTime_long, &maxTime_ensemble, &intervallini, &N_series, &n_bins, &x0,
                            &seed_ensemble, &seed_long, &max_lag_seconds_from_file, &ensemble_autocorr_toggle);
    fclose(infile);
    //l'autocorrelation toggle spegne e accende il calcolo dell'autocorrelazione per l'insieme.
    //Quando è acceso (impostare 1 nell'input) il programma potrebbe richiedere molto tempo'

    if (numeri_input < 10) {
        fprintf(stderr, "Err: 'input_secondo.dat' deve contenere 10 parametri.\n");
        fprintf(stderr, "Ricorda: interruttore autocorr dell'ensemble alla fine (0=attiva, 1=disattiva).\n");
        return 1;
    }

    if (maxTime_ensemble > 0 && intervallini > 0 && maxTime_ensemble > INT_MAX / intervallini) {
        fprintf(stderr, "Attenzione: I valori di input per la media d'insieme causano un integer overflow.\n");
        return 1;
    }

    if (max_lag_seconds_from_file < 0) {
        fprintf(stderr, "Attenzione: max_lag_seconds dall'input è negativo. Imposto a 0.\n");
        max_lag_seconds_from_file = 0;
    }

    // calcolo separato del numero degli step per serie lunga e serie piccole
    long long int long_N_steps = (long long int) maxTime_long * intervallini;
    long long int short_N_steps = (long long int)maxTime_ensemble * intervallini;
    if (short_N_steps == 0) short_N_steps = 1;
    double DT = 1.0 / intervallini;

    // Controllo dei parametri e vado ad allocare la memoria che mi serve
    printf("Inizio: Media d'insieme\n");
    printf("Generazione di %d serie brevi da %lld punti ciascuna.\n", N_series, short_N_steps);
    printf("Controllo parametri simulazione insieme:\n");
    printf("Tempo fisico: %d\n", maxTime_ensemble);
    printf("N serie: %d\n", N_series);
    printf("N bin ist.: %d\n", n_bins);
    printf("seed input: %d\n", seed_ensemble);
    printf("x0 e' generato gaussianamente per ogni serie.\n\n");


    if (seed_ensemble <= 0) srand((unsigned int)time(NULL));
    else srand((unsigned int)seed_ensemble);


    //Apro file per la stampa delle serie
    FILE *ensemble_paths_file = fopen("ensemble_paths.dat", "w");
    if (ensemble_paths_file == NULL) {
        perror("Errore nell'apertura di ensemble_paths.dat");
        return 1;
    }

    size_t global_ensemble_size = (size_t)N_series * short_N_steps;
    struct Vector *global_ensemble_data = newVector(global_ensemble_size);
    if(global_ensemble_data == NULL) {
        perror("Allocazione memoria fallita per global_ensemble_data");
        fclose(ensemble_paths_file);
        return 1;
    }
    size_t data_cont = 0;

    struct Vector **all_short_series = (struct Vector **)malloc(N_series * sizeof(struct Vector *));
    if(all_short_series == NULL) {
        perror("Allocazione memoria fallita per all_short_series");
        delVector(global_ensemble_data);
        fclose(ensemble_paths_file);
        return 1;
    }

    for (int i = 0; i < N_series; i++) {
        all_short_series[i] = newVector(short_N_steps);
        if(all_short_series[i] == NULL) {
            perror("Allocazione memoria fallita generazione all_short_series");
            for(int k=0; k<i; ++k) delVector(all_short_series[k]);
            free(all_short_series);
            delVector(global_ensemble_data);
            fclose(ensemble_paths_file);
            return 1;
        }

        /*
        Ciclo for di integrazione della SDE.
        Genero gaussianamente un punto iniziale per le n serie che sto integrando.
        if-else in forma short-hand dato che il drift ha una forma semplice.
        Uso lo schema di Eulero per integrare la traiettoria, il rumore è generato gaussiano
        e scalato tramite la radice del passo d'integrazione.
        */


        all_short_series[i]->data[0] = box_muller();

        for (long long int j = 1; j < short_N_steps; j++) {
            double x = all_short_series[i]->data[j-1];
            double h_x = (x >= 0) ? -K : K; //
            double noise = sqrt(2.0 * D * DT) * box_muller();
            all_short_series[i]->data[j] = x + h_x * DT + noise;
        }
    }


    /*Trasferisco tutte le serie generate in un file globale, che viene campionato ogni secondo reale
     * e stampato su file per poter eventualmente controllare il profilo delle traiettorie generate.
     */

    for (int i = 0; i < N_series; i++)
    {
        for (long long int j = 0; j < short_N_steps; j++)
        {
            global_ensemble_data->data[data_cont++] = all_short_series[i]->data[j];

            if (j % intervallini == 0)
            {
                fprintf(ensemble_paths_file, "%.6f ", all_short_series[i]->data[j]);
            }
        }
        fprintf(ensemble_paths_file, "\n");
    }
    fclose(ensemble_paths_file);

    //Gestione dell'analisi delle traiettorie per le medie d'insieme

    if (ensemble_autocorr_toggle == 1) {
        printf("Inizio: Calcolo autocorrelazione per l'ensemble\n");
        calcolo_autocorr_ens_average(all_short_series, N_series, maxTime_ensemble, intervallini, "autocorrelation_ensemble.dat");
        printf("Fine: Autocorrelazione dell'ensemble stampata\n");
    } else {
        printf("Calcolo autocorrelazione per l'ensemble spento (interruttore = %d).\n", ensemble_autocorr_toggle);
    }

    //Libero la memoria del vettore che contiene anche tutte le frazioni di secondo
    for(int i = 0; i < N_series; i++){
        delVector(all_short_series[i]);
    }
    free(all_short_series);

    //Chiamo tutte le funzioni di analisi
    printf("Inizio analisi ensemble\n");
    FILE *moments_ensemble_file = fopen("moments_ensemble.dat", "w");
    if (moments_ensemble_file == NULL) {
        perror("Errore nell'apertura di moments_ensemble.dat");
        delVector(global_ensemble_data);
        return 1;
    }
    istogramma(global_ensemble_data, n_bins, "ensemble_globale");
    calcolo_stampa_momenti(global_ensemble_data, "serie piccole (ensemble)", moments_ensemble_file);
    fclose(moments_ensemble_file);
    delVector(global_ensemble_data);
    printf("Fine analisi medie d'insieme.\n\n");


    //Calcoli per la serie temporale lunga
    //sposto sotto il main per leggibilità

    printf("Inizio: Serie Temporale Lunga\n");
    printf("Controllo parametri simulazione lunga:\n");
    printf("Tempo fisico: %d\n", maxTime_long);
    printf("Intervalli per unita' di tempo: %d\n", intervallini);
    printf("Numero totale di step: %lld\n", long_N_steps);
    printf("x0 input: %.2f\n", x0);
    printf("seed: %d\n", seed_long);
    printf("Massimo lag da input per autocorrelazione: %d\n\n", max_lag_seconds_from_file);

    run_long_timeseries_analysis(long_N_steps, DT, x0, seed_long, n_bins, max_lag_seconds_from_file);
    printf("Fine analisi serie lunga.\n\n");

    printf("Fine scrittura file output.\n");
    return 0;
}


//SERIE TEMP LUNGA
//passo gli input di file che mi servono per avviare la simulazione
void run_long_timeseries_analysis(long long int n_steps, double dt, double x0, int seed, int n_bins, int max_lag_seconds) {

    const int DOWNSAMPLE_FACTOR = (int) (1.0 / dt); //recupero il numero intervallini

    //calcolo quanti punti della simulazione (secondi) vado ad immagazzinare
    //alleggerisco  l'uso della memoria salvando solo i secondi
    //teoricamente posso arrivare a serie temporali che simulano "decine di anni" senza eccedere 16gb di ram
    //ma a quel punto il programma richiede ~1hr per terminare la simulazione
    size_t n_points_to_store = n_steps / DOWNSAMPLE_FACTOR;
    if (n_steps > 0 && n_steps % DOWNSAMPLE_FACTOR != 0) {
        n_points_to_store++;
    }

    if (n_points_to_store == 0 && n_steps > 0) {
        n_points_to_store = 1;
    }
    if (n_steps == 0) {
        printf("Numero di steps per la serie lunga è zero. Nessuna simulazione.\n");
        return;
    }


    printf("Simulazione di %lld passi. Salvo in memoria %zu punti (downsample %d).\n",
           n_steps, n_points_to_store, DOWNSAMPLE_FACTOR);

    if (seed <= 0) srand((unsigned int)time(NULL) + 1);
    else srand((unsigned int)seed);

    struct Vector* long_series_downsampled = newVector(n_points_to_store);
    if (long_series_downsampled == NULL) {
        perror("Allocazione di memoria fallita per la serie lunga campionata");
        return;
    }

    double x_current = x0;
    size_t cont = 0;

    //Integrazione della traiettoria
    for (long long int j = 0; j < n_steps; j++) {
        if (j % DOWNSAMPLE_FACTOR == 0 || j == n_steps - 1) {
            if (cont < n_points_to_store) {
                long_series_downsampled->data[cont] = x_current;
                cont++;
            }
        }

        double h_x = (x_current >= 0) ? -K : K;
        double noise = sqrt(2.0 * D * dt) * box_muller();
        x_current = x_current + h_x * dt + noise;
    }

    printf("Inizio: analisi serie temporale lunga\n");

    FILE *moments_long_file = fopen("moments_long_series.dat", "w");
    if (moments_long_file == NULL) {
        perror("Errore nell'apertura di moments_long_series.dat");
        delVector(long_series_downsampled);
        return;
    }
    istogramma(long_series_downsampled, n_bins, "long_series_globale");
    calcolo_stampa_momenti(long_series_downsampled, "Serie temp (campionata)", moments_long_file);
    fclose(moments_long_file);

    FILE *long_ts_file = fopen("timeseries_path.dat", "w");
    if(long_ts_file != NULL){
        for(size_t i = 0; i < long_series_downsampled->size; i++){
            fprintf(long_ts_file, "%.6f\n", long_series_downsampled->data[i]);
        }
        fclose(long_ts_file);
    } else {
        perror("Errore nell'apertura di timeseries_path.dat");
    }


    printf("Calcolo autocorrelazione serie lunga\n");
    calcolo_autocorr_ts(long_series_downsampled, max_lag_seconds, "autocorrelation_long_series.dat");
    printf("Autocorrelazione stampata su file. FINE\n");

    delVector(long_series_downsampled);
}


// funzioni ausiliarie

//Generazione numeri random generati gaussianamente

double box_muller() {
    static double z1;
    static int generate = 0;
    generate = !generate;
    if (!generate) return z1;
    double u1, u2;
    do {
        u1 = rand() * (1.0 / RAND_MAX);
        u2 = rand() * (1.0 / RAND_MAX);
    } while (u1 <= 1e-9);
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    z1 = sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);
    return z0;
}

//Colloco i punti generati nelle traiettorie in un istogramma diviso linearmente

void istogramma (struct Vector * array, int binNum, const char *name) {
    if (array == NULL || array->size == 0) {
        fprintf(stderr, "Impossibile creare istogramma: vettore nullo o vuoto per %s.\n", name);
        return;
    }
    char * filename = malloc(sizeof(char)*50);
    if (filename == NULL) {
        perror("Allocazione memoria fallita per nome file istogramma");
        return;
    }
    snprintf(filename, 50, "istogramma_%s.dat", name);
    double massimo = max(array);
    double minimo = min(array);

    if (fabs(massimo - minimo) < 1e-9) {
        fprintf(stderr, "Tutti i valori sono identici. Istogramma non generato per %s.\n", name);
        free(filename);
        return;
    }

    double binSize = (massimo - minimo) / binNum;


    struct Vector * istogramma_counts = newVector(binNum);
    if (istogramma_counts == NULL) {
        perror("Allocazione memoria fallita per istogramma_counts");
        free(filename);
        return;
    }
    for(size_t i = 0; i < istogramma_counts->size; ++i) istogramma_counts->data[i] = 0;


    for(size_t i=0; i < array->size; i++) {
        int bin = (int)((array->data[i] - minimo) / binSize);
        if (bin >= binNum) bin = binNum - 1;
        if (bin < 0) bin = 0;
        istogramma_counts->data[bin]++;
    }
    long long total_counts = 0;
    for(int i=0; i<binNum; i++){
        total_counts += (long long)istogramma_counts->data[i];
    }
    printf("Controllo Istogramma: Conteggi totali = %lld (su %zu dati) \n", total_counts, array->size);

    FILE *fp4 = fopen(filename, "w");
    if (fp4 == NULL) {
        perror("Errore nell'apertura del file istogramma");
        free(filename);
        delVector(istogramma_counts);
        return;
    }

    fprintf(fp4, "# Centro_Bin\tDensita_di_Probabilita\n");
    double total_normalized_area = 0.0;
    for(int i=0; i < binNum; i++) {
        double bin_center = minimo + (i + 0.5) * binSize;
        double normalized_height = (binSize > 0) ? (istogramma_counts->data[i] / (double)array->size) / binSize : 0.0;
        total_normalized_area += normalized_height * binSize;
        fprintf(fp4, "%f\t%f\n", bin_center, normalized_height);
    }

    printf("Controllo Istogramma 2: L'area totale calcolata e' %f.\n", total_normalized_area);
    fclose(fp4);
    free(filename);
    delVector(istogramma_counts);
}


//Calcolo 1o 2o 3o 4o momento e varianza della serie

void calcolo_stampa_momenti(struct Vector* vec, const char* series_type, FILE *outfile) {
    if (vec == NULL || vec->size == 0 || outfile == NULL) {
        fprintf(stderr, "Impossibile calcolare i momenti: vettore nullo o vuoto o file non valido per %s.\n", series_type);
        return;
    }
    double m1 = 0.0, m2 = 0.0, m3 = 0.0, m4 = 0.0;
    for (size_t i = 0; i < vec->size; i++) {
        double x = vec->data[i];
        m1 += x;
        m2 += x * x;
        m3 += x * x * x;
        m4 += x * x * x * x;
    }
    m1 /= vec->size;
    m2 /= vec->size;
    m3 /= vec->size;
    m4 /= vec->size;
    double varianza = m2 - (m1 * m1);
    fprintf(outfile, "Momenti %s\n", series_type);
    fprintf(outfile, "Media <x>:\t\t%.6f\n", m1);
    fprintf(outfile, "Media Quadratica <x^2>:\t%.6f\n", m2);
    fprintf(outfile, "Momento Terzo <x^3>:\t%.6f\n", m3);
    fprintf(outfile, "Momento Quarto <x^4>:\t%.6f\n", m4);
    fprintf(outfile, "Varianza (<x^2>-<x>^2):\t%.6f\n\n", varianza);
}


double max (struct Vector * array) {
    if (array == NULL || array->size == 0) {
        fprintf(stderr, "Errore: input nullo o vuoto per la funzione max.\n");
        return 0.0;
    }
    double maxTemp = array->data[0];
    for (size_t i=1; i < array->size; i++) {
        if (array->data[i] > maxTemp) {
            maxTemp = array->data[i];
        }
    }
    return maxTemp;
}

double min (struct Vector * array) {
    if (array == NULL || array->size == 0) {
        fprintf(stderr, "Errore: input nullo o vuoto per la funzione min.\n");
        return 0.0;
    }
    double minTemp = array->data[0];
    for (size_t i=1; i < array->size; i++) {
        if (array->data[i] < minTemp) {
            minTemp = array->data[i];
        }
    }
    return minTemp;
}


//Calcolo autocorrelazione dedicato per le serie di ensemble
//Al massimo mi concentro sull'1% della serie generata, dopo qualche simulazione ho imposto un hard limit al lag dell'autocorrelazione
//sommo la covarianza ad un certo lag, la normalizzo, divido per il numero di serie e memorizzo in un vettore con tutti i lag

void calcolo_autocorr_ens_average(struct Vector** all_series, int n_series, int maxTime_ensemble, int intervallini, const char* filename) {
    if (all_series == NULL || n_series == 0 || intervallini <= 0) {
        fprintf(stderr, "Input non valido per l'autocorrelazione dell'ensemble.\n");
        return;
    }

    // New logic: calculate lag as time/100, but cap it at 50.
    int max_lag_seconds = maxTime_ensemble / 100;
    if (max_lag_seconds > 30) {
        max_lag_seconds = 30;
    }
    if (max_lag_seconds < 1 && maxTime_ensemble > 0) {
        max_lag_seconds = 1;
    }

    printf("Autocorrelazione ens: uso un lag massimo di %d secondi.\n", max_lag_seconds);

    struct Vector* autocorr_sum = newVector(max_lag_seconds + 1);
    if (autocorr_sum == NULL) {
        perror("Allocazione memoria fallita per autocorr_sum");
        return;
    }
    for (int i = 0; i <= max_lag_seconds; ++i) {
        autocorr_sum->data[i] = 0.0;
    }

    int series_num = 0;

    // ciclo for su tutte le serie dell'insieme'
    for (int i = 0; i < n_series; i++) {
        printf("\rElaborazione serie ensemble %d/%d", i + 1, n_series);
        fflush(stdout);

        struct Vector* current_series = all_series[i];
        if (current_series == NULL || current_series->size == 0) continue;

        double media = 0.0;
        for (size_t j = 0; j < current_series->size; j++) media += current_series->data[j];
        media /= current_series->size;

        double varianza = 0.0;
        for (size_t j = 0; j < current_series->size; j++) varianza += pow(current_series->data[j] - media, 2);
        varianza /= current_series->size;

        if (fabs(varianza) < 1e-12) continue;

        series_num++;

        // ciclo for sui secondi del lag
        for (int lag_s = 0; lag_s <= max_lag_seconds; lag_s++) {
            // trasformo l'input per puntarlo sui secondi pieni'
            long long int tau = (long long int)lag_s * intervallini;

            if (tau >= (long long int)current_series->size) {
                break;
            }

            double cov = 0.0;
            long long int num_coppie = current_series->size - tau;

            for (size_t k = 0; k < (size_t)num_coppie; k++) {
                cov += (current_series->data[k] - media) * (current_series->data[k + tau] - media);
            }
            cov /= num_coppie;

            double normalized_corr = cov / varianza;
            autocorr_sum->data[lag_s] += normalized_corr;
        }
    }

    printf("\n");

    if (series_num == 0) {
        printf("Nessuna serie valida per il calcolo dell'autocorrelazione dell'ensemble.\n");
        delVector(autocorr_sum);
        return;
    }

    FILE* fp = fopen(filename, "w");
    if (fp == NULL) {
        perror("Errore nell'apertura del file per l'autocorrelazione dell'ensemble");
        delVector(autocorr_sum);
        return;
    }

    fprintf(fp, "# Lag_Time(s)\tAverage_Autocorrelation\n");

    /////scrittura
    for (int lag_s = 0; lag_s <= max_lag_seconds; lag_s++) {
        double avg_corr = autocorr_sum->data[lag_s] / series_num;
        double lag_ind = (double)lag_s;
        fprintf(fp, "%.6f\t%.6f\n", lag_ind, avg_corr);
    }

    fclose(fp);
    delVector(autocorr_sum);
}


//Funzione autocorrelazione che uso per la serie lunga

void calcolo_autocorr_ts(struct Vector* vec, int max_lag_seconds, const char* filename) {
    if (vec == NULL || vec->size == 0) {
        printf("Impossibile calcolare l'autocorrelazione, vettore non valido.\n");
        return;
    }

    double media = 0.0;
    for (size_t i = 0; i < vec->size; i++) {
        media += vec->data[i];
    }
    media /= vec->size;

    double varianza = 0.0;
    for (size_t i = 0; i < vec->size; i++) {
        varianza += (vec->data[i] - media) * (vec->data[i] - media);
    }
    varianza /= vec->size;

    if (fabs(varianza) < 1e-12) {
        printf("Varianza nulla, autocorrelazione non definita.\n");
        return;
    }

    long long int max_lag_indici_campionato = max_lag_seconds;

    if (max_lag_indici_campionato >= (long long int)vec->size) {
        max_lag_indici_campionato = (long long int)vec->size - 1;
        printf("ATTENZIONE: Massimo lag (%d secondi) supera la lunghezza della serie campionata.\n",
               max_lag_seconds, vec->size);
        printf("Il calcolo dell'autocorrelazione sarà limitato a %lld secondi.\n", max_lag_indici_campionato);
    }
    if (max_lag_indici_campionato < 0) max_lag_indici_campionato = 0;


    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        perror("Errore nell'apertura del file per l'autocorrelazione");
        return;
    }

    fprintf(fp, "# Lag_Time(s)\tAutocorrelation\n");

    for (long long int tau = 0; tau <= max_lag_indici_campionato; tau++) {
        double covarianza_at_lag = 0.0;
        long long int num_coppie_for_lag = vec->size - tau;

        if (num_coppie_for_lag <= 0) {
            break;
        }

        for (size_t i = 0; i < (size_t)num_coppie_for_lag; i++) {
            covarianza_at_lag += (vec->data[i] - media) * (vec->data[i + tau] - media);
        }
        covarianza_at_lag /= num_coppie_for_lag;

        double autocorrelation_coefficient = covarianza_at_lag / varianza;

        double lag_ind = (double)tau;

        fprintf(fp, "%.6f\t%.6f\n", lag_ind, autocorrelation_coefficient);
    }

    fclose(fp);
}
