#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
//
// DOVE NON COMMENTATO VEDERE CODICE K-DRIFT
// qui cerco di commentare solo quanto è diverso dal primo programma
//
// costanti
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
    if (retVal->data == NULL) {
        free(retVal);
        return NULL;
    }
    retVal->size = sz;
    return retVal;
}

void delVector(struct Vector *vector) {
    if (vector != NULL) {
        free(vector->data);
        free(vector);
    }
}

// funzioni prototipo
double box_muller();
void calcolo_stampa_momenti(struct Vector* vec, const char* series_type, FILE *outfile);
void calcolo_autocorr_ts(struct Vector* vec, const char* filename, int max_lag);
double max(struct Vector * array);
double min(struct Vector * array);
void istogramma_lineare(struct Vector * array, int binNum, const char *name);
void istogramma_logaritmico(struct Vector * array, int binNum, const char *name);
void run_long_timeseries_analysis_alpha(long long int n_steps, double dt, double x0, double alpha, int seed, int n_bins, int downsample_factor, int max_lag);
void calcolo_autocorr_ens_average(struct Vector** all_series, int n_series, int maxTime_ensemble, int intervallini, const char* filename);

//
// main
//

int main() {
    // dichiarazioni e lettura da file
    int maxTime_long, maxTime_ensemble, intervallini, N_series;
    int seed_ensemble, n_bins, seed_long, max_lag, ensemble_autocorr_toggle;
    double x0, alpha;

    FILE *infile = fopen("input_secondo_alpha.dat", "r");
    if (infile == NULL) {
        perror("Err: non esiste 'input_secondo_alpha.dat'");
        return 1;
    }

    int numeri_input = fscanf(infile, "%d %d %d %d %d %lf %lf %d %d %d %d",
                            &maxTime_long, &maxTime_ensemble, &intervallini,
                            &N_series, &n_bins, &x0, &alpha,
                            &seed_ensemble, &seed_long, &max_lag, &ensemble_autocorr_toggle);
    fclose(infile);

    if (numeri_input < 11) {
        fprintf(stderr, "Err: 'input_secondo_alpha.dat' non contiene abbastanza valori. Sono richiesti 11 parametri.\n");
        fprintf(stderr, "Assicurati di includere il toggle per l'autocorrelazione dell'ensemble alla fine (1=attivata, 0=disattivata).\n");
        return 1;
    }

    if (max_lag <= 0) {
        fprintf(stderr, "Attenzione: 'max_lag' deve essere un intero positivo.\n");
        return 1;
    }

    if (fabs(x0) < 1e-9) {
        fprintf(stderr, "Attenzione: x0 e' molto vicino a zero. Scegliere un x0 diverso da zero.\n");
    }

    if (intervallini <= 0) {
        fprintf(stderr, "Errore: 'intervallini' deve essere un intero positivo.\n");
        return 1;
    }

    if (maxTime_ensemble > 0 && intervallini > 0 && maxTime_ensemble > INT_MAX / intervallini) {
        fprintf(stderr, "Attenzione: I valori di input per la media d'insieme causano un integer overflow.\n");
        return 1;
    }

    long long int long_N_steps = (long long int) maxTime_long * intervallini;
    long long int short_N_steps = (long long int) maxTime_ensemble * intervallini;
    if (short_N_steps == 0) short_N_steps = 1;
    double DT = 1.0 / intervallini;

    // Calcolo Ensemble
    printf("Media d'insieme (alpha drift)\n");
    printf("Generazione di %d serie brevi da %lld punti ciascuna.\n", N_series, short_N_steps);
    printf("x0 = Generato gaussianamente per ogni serie.\n\n");


    if (seed_ensemble <= 0) srand((unsigned int)time(NULL));
    else srand((unsigned int)seed_ensemble);

    size_t global_ensemble_size = (size_t)N_series * short_N_steps;
    struct Vector *global_ensemble_data = newVector(global_ensemble_size);
    if(global_ensemble_data == NULL) { perror("Allocazione memoria fallita"); return 1; }
    size_t data_cont = 0;

    struct Vector **all_short_series = (struct Vector **)malloc(N_series * sizeof(struct Vector *));
    if(all_short_series == NULL) { perror("Allocazione memoria fallita"); return 1; }

    /*
    In questo caso il drift si accende e spegne oltre una certa distanza dall'origine
    per evitare di avere problemi nella divisione all'interno del calcolo del drift
    nel codice a-drift regolarizzato provo ad affrontare il problema cercando un drift continuo
    ma con comportamento asintotico simile a questo
    */

    for (int i = 0; i < N_series; i++) {
        all_short_series[i] = newVector(short_N_steps);
        if(all_short_series[i] == NULL) { perror("Allocazione memoria fallita"); return 1; }

        all_short_series[i]->data[0] = box_muller();

        for (long long int j = 1; j < short_N_steps; j++) {
            double x = all_short_series[i]->data[j-1];
            double h_x = 0.0;
            if (fabs(x) > 1e-9) {
                h_x = -alpha / x;
            }
            double noise = sqrt(2.0 * D * DT) * box_muller();
            all_short_series[i]->data[j] = x + h_x * DT + noise;
        }
    }

    FILE *ensemble_paths_file = fopen("ensemble_paths_alpha.dat", "w");
    for (int i = 0; i < N_series; i++) {
        for (long long int j = 0; j < short_N_steps; j++) {
            global_ensemble_data->data[data_cont++] = all_short_series[i]->data[j];
            if (j % intervallini == 0 && ensemble_paths_file != NULL) {
                fprintf(ensemble_paths_file, "%.6f ", all_short_series[i]->data[j]);
            }
        }
        if (ensemble_paths_file != NULL) fprintf(ensemble_paths_file, "\n");
    }
    if (ensemble_paths_file != NULL) fclose(ensemble_paths_file);

    if (ensemble_autocorr_toggle == 1) {
        printf("Calcolo autocorrelazione per l'ensemble abilitato\n");
        calcolo_autocorr_ens_average(all_short_series, N_series, maxTime_ensemble, intervallini, "autocorrelation_ensemble_alpha.dat");
        printf("Autocorrelazione dell'ensemble stampata.\n");
    } else {
        printf("Calcolo autocorrelazione per l'ensemble spento (toggle = %d).\n", ensemble_autocorr_toggle);
    }


    for(int i = 0; i < N_series; i++){
        delVector(all_short_series[i]);
    }
    free(all_short_series);

    printf("Inizio analisi ensemble\n");
    FILE *moments_ensemble_file = fopen("moments_ensemble_alpha.dat", "w");
    istogramma_lineare(global_ensemble_data, n_bins, "ensemble_globale_alpha");
    istogramma_logaritmico(global_ensemble_data, n_bins, "ensemble_globale_alpha");
    calcolo_stampa_momenti(global_ensemble_data, "serie piccole (ensemble Alpha)", moments_ensemble_file);
    fclose(moments_ensemble_file);
    delVector(global_ensemble_data);
    printf("Fine analisi medie d'insieme.\n\n");

    // Calcoli per la serie temporale lunga, vedi sotto
    printf("Serie Temporale Lunga (alpha drift)\n");
    printf("x0 fisso: %.2f\n", x0);
    run_long_timeseries_analysis_alpha(long_N_steps, DT, x0, alpha, seed_long, n_bins, intervallini, max_lag);
    printf("Fine analisi serie lunga.\n\n");

    printf("Fine scrittura file output.\n");
    return 0;
}


// SERIE TEMP LUNGA
void run_long_timeseries_analysis_alpha(long long int n_steps, double dt, double x0, double alpha, int seed, int n_bins, int downsample_factor, int max_lag) {

    size_t n_points_to_store = n_steps / downsample_factor;
    if (n_points_to_store == 0 && n_steps > 0) n_points_to_store = 1;
    if (n_steps == 0) return;

    printf("Simulazione di %lld passi. Salvo in memoria %zu punti (1 ogni %d passi).\n",
           n_steps, n_points_to_store, downsample_factor);

    if (seed <= 0) srand((unsigned int)time(NULL) + 1);
    else srand((unsigned int)seed);

    struct Vector* long_series_downsampled = newVector(n_points_to_store);
    if (long_series_downsampled == NULL) {
        perror("Allocazione memoria fallita per la serie lunga campionata");
        return;
    }

    double x_current = x0;
    size_t cont = 0;

    for (long long int j = 0; j < n_steps; j++) {
        if (j % downsample_factor == 0) {
            if (cont < n_points_to_store) {
                long_series_downsampled->data[cont] = x_current;
                cont++;
            }
        }

        double h_x = 0.0;
        if (fabs(x_current) > 1e-9) {
            h_x = -alpha / x_current;
        }

        double noise = sqrt(2.0 * D * dt) * box_muller();
        x_current = x_current + h_x * dt + noise;
    }

    printf("Inizio analisi sulla serie temporale lunga\n");
    FILE *moments_long_file = fopen("moments_long_series_alpha.dat", "w");

    istogramma_lineare(long_series_downsampled, n_bins, "long_series_globale_alpha");
    istogramma_logaritmico(long_series_downsampled, n_bins, "long_series_globale_alpha");
    calcolo_stampa_momenti(long_series_downsampled, "Serie temporale (campionata, alpha)", moments_long_file);

    printf("Inizio calcolo autocorrelazione TS\n");
    calcolo_autocorr_ts(long_series_downsampled, "autocorrelation_alpha.dat", max_lag);
    printf("Fine calcolo autocorrelazione.\n");

    fclose(moments_long_file);

    FILE *long_ts_file = fopen("timeseries_path_alpha.dat", "w");
    if(long_ts_file != NULL){
        for(size_t i = 0; i < long_series_downsampled->size; i++){
            fprintf(long_ts_file, "%.6f\n", long_series_downsampled->data[i]);
        }
        fclose(long_ts_file);
    }

    delVector(long_series_downsampled);
}


// Funzioni ausiliarie


//In questo caso il limite hard-coded per l'autocorrelazione è 200, dopo questo il programma diventa troppo lento.
//Dopo diverse prove 200 penso sia lo sweet-spot tra l'avere un'idea della memoria del processo e la velocità computazionale
void calcolo_autocorr_ens_average(struct Vector** all_series, int n_series, int maxTime_ensemble, int intervallini, const char* filename) {
    if (all_series == NULL || n_series == 0 || intervallini <= 0) {
        fprintf(stderr, "Input non valido per l'autocorrelazione dell'ensemble.\n");
        return;
    }

    int max_lag_seconds = maxTime_ensemble / 100;
    if (max_lag_seconds > 200) {
        max_lag_seconds = 200;
    }
    if (max_lag_seconds < 1 && maxTime_ensemble > 0) {
        max_lag_seconds = 1;
    }

    printf("Autocorrelazione ensemble: uso un lag massimo di %d secondi.\n", max_lag_seconds);

    struct Vector* autocorr_sum = newVector(max_lag_seconds + 1);
    if (autocorr_sum == NULL) {
        perror("Allocazione memoria fallita per autocorr_sum");
        return;
    }
    for (int i = 0; i <= max_lag_seconds; ++i) {
        autocorr_sum->data[i] = 0.0;
    }

    int series_num = 0;

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

        for (int lag_s = 0; lag_s <= max_lag_seconds; lag_s++) {
            long long int tau = (long long int)lag_s * intervallini;
            if (tau >= (long long int)current_series->size) break;

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
        printf("Nessuna serie valida trovata per il calcolo dell'autocorrelazione dell'ensemble.\n");
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

    for (int lag_s = 0; lag_s <= max_lag_seconds; lag_s++) {
        double avg_corr = autocorr_sum->data[lag_s] / series_num;
        double lag_ind = (double)lag_s;
        fprintf(fp, "%.6f\t%.6f\n", lag_ind, avg_corr);
    }

    fclose(fp);
    delVector(autocorr_sum);
}


void calcolo_autocorr_ts(struct Vector* vec, const char* filename, int max_lag) {
    if (vec == NULL || vec->size < 2) {
        fprintf(stderr, "Impossibile calcolare l'autocorrelazione, vettore non valido\n");
        return;
    }

    FILE* outfile = fopen(filename, "w");
    if (outfile == NULL) {
        perror("Errore nell'apertura del file per l'autocorrelazione");
        return;
    }

    if (max_lag >= (int)vec->size) {
        printf("Attenzione: max_lag (%d) e' maggiore o uguale della dimensione della serie (%zu). Verra' usato un max_lag di %zu.\n", max_lag, vec->size, vec->size - 1);
        max_lag = vec->size - 1;
    }

    printf("Calcolo autocorrelazione per %d lag.\n", max_lag);

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
        fprintf(stderr, "Varianza nulla, autocorrelazione non definita.\n");
        fclose(outfile);
        return;
    }

    fprintf(outfile, "# Lag\tAutocorrelation_Normalized\n");

    for (int lag = 0; lag < max_lag; lag++) {
        double C_lag = 0.0;
        for (size_t i = 0; i < vec->size - lag; i++) {
            C_lag += (vec->data[i] - media) * (vec->data[i + lag] - media);
        }
        C_lag /= (vec->size - lag);

        double normalized_C = C_lag / varianza;
        fprintf(outfile, "%d\t%.6f\n", lag, normalized_C);
    }

    fclose(outfile);
}


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
/*
In questo codice abbiamo la possibilità di usare due tipi di istogramma
Istogramma lineare e istogramma logaritmico
Ho scelto di utilizzare il binning logaritmico a causa della grande
variabilità del processo, che è molto soggetto agli eventi estremi
a causa della forma del drift. Con un istogramma di tipo logaritmico
si mitiga molto la rumorosità delle code, che contengono un numero di dati
molto piccolo nel caso di un binning lineare.
*/


void istogramma_lineare(struct Vector * array, int binNum, const char *name) {
    if (array == NULL || array->size == 0) return;
    char * filename = malloc(sizeof(char)*100);
    if (filename == NULL) return;
    snprintf(filename, 100, "istogramma_lineare_%s.dat", name);

    double massimo = max(array);
    double minimo = min(array);
    double binSize = (massimo - minimo) / binNum;

    struct Vector * istogramma_counts = newVector(binNum);
    if (istogramma_counts == NULL) { free(filename); return; }
    for(size_t i = 0; i < istogramma_counts->size; ++i) istogramma_counts->data[i] = 0;

    for(size_t i=0; i < array->size; i++) {
        int bin = (int)((array->data[i] - minimo) / binSize);
        if (bin >= binNum) bin = binNum - 1;
        if (bin < 0) bin = 0;
        istogramma_counts->data[bin]++;
    }

    printf("Controllo Istogramma Lineare: ...\n");
    FILE *fp4 = fopen(filename, "w");
    fprintf(fp4, "# Centro_Bin\tDensita_di_Probabilita\n");
    for(int i=0; i < binNum; i++) {
        double bin_center = minimo + (i + 0.5) * binSize;
        double normalized_height = (binSize > 0) ? (istogramma_counts->data[i] / (double)array->size) / binSize : 0;
        fprintf(fp4, "%f\t%f\n", bin_center, normalized_height);
    }
    fclose(fp4);
    free(filename);
    delVector(istogramma_counts);
}


//Istogramma che ho usato nelle simulazioni "ufficiali":
void istogramma_logaritmico(struct Vector * array, int binNum, const char *name) {
    if (array == NULL || array->size == 0) return;

    double min_abs = DBL_MAX;
    double max_abs = 0.0;
    for (size_t i = 0; i < array->size; i++) {
        double val_abs = fabs(array->data[i]);
        if (val_abs > 0) {
            if (val_abs < min_abs) min_abs = val_abs;
            if (val_abs > max_abs) max_abs = val_abs;
        }
    }

    if (min_abs >= max_abs) {
        printf("Impossibile creare un istogramma logaritmico: range di dati non valido.\n");
        return;
    }

    struct Vector* istogramma_counts = newVector(binNum);
    if (istogramma_counts == NULL) return;
    for (size_t i = 0; i < istogramma_counts->size; ++i) istogramma_counts->data[i] = 0;

    double log_min = log(min_abs);
    double log_max = log(max_abs);
    double log_bin_size = (log_max - log_min) / binNum;

    size_t punti_totali_istogramma = 0;
    for (size_t i = 0; i < array->size; i++) {
        double val_abs = fabs(array->data[i]);
        if (val_abs > 0) {
            int bin_index = (int)((log(val_abs) - log_min) / log_bin_size);
            if (bin_index >= 0 && bin_index < binNum) {
                istogramma_counts->data[bin_index]++;
                punti_totali_istogramma++;
            }
        }
    }

    char nome_file[100];
    snprintf(nome_file, 100, "istogramma_logaritmico_%s.dat", name);
    FILE *file = fopen(nome_file, "w");
    if (file == NULL) {
        delVector(istogramma_counts);
        return;
    }

    fprintf(file, "# CentroBin\tDensitaDiProbabilita\n");

    for (int i = 0; i < binNum; i++) {
        double conteggio = istogramma_counts->data[i];

        if (conteggio > 0) {
            double log_limite_inf = log_min + i * log_bin_size;
            double log_limite_sup = log_min + (i + 1.0) * log_bin_size;

            double limite_inf_lin = exp(log_limite_inf);
            double limite_sup_lin = exp(log_limite_sup);

            double ampiezza_bin_lin = limite_sup_lin - limite_inf_lin;
            double centro_bin = sqrt(limite_inf_lin * limite_sup_lin);

            double valore_normalizzato = 0.0;
            if (punti_totali_istogramma > 0 && ampiezza_bin_lin > 0) {
                valore_normalizzato = conteggio / ((double)punti_totali_istogramma * ampiezza_bin_lin);
            }

            fprintf(file, "%e\t%e\n", centro_bin, valore_normalizzato);
        }
    }

    fclose(file);
    delVector(istogramma_counts);
    printf("Istogramma normalizzato salvato in %s\n", nome_file);
}


void calcolo_stampa_momenti(struct Vector* vec, const char* series_type, FILE *outfile) {
    if (vec == NULL || vec->size == 0 || outfile == NULL) return;
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
    if (array == NULL || array->size == 0) return 0.0;
    double maxTemp = array->data[0];
    for (size_t i=1; i < array->size; i++) {
        if (array->data[i] > maxTemp) {
            maxTemp = array->data[i];
        }
    }
    return maxTemp;
}

double min (struct Vector * array) {
    if (array == NULL || array->size == 0) return 0.0;
    double minTemp = array->data[0];
    for (size_t i=1; i < array->size; i++) {
        if (array->data[i] < minTemp) {
            minTemp = array->data[i];
        }
    }
    return minTemp;
}
