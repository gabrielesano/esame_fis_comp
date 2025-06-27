#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <limits.h>

//
//Codice di controllo per la prova di concetto della parallelizzazione
//Esegue le stesse operazioni del codice parallelo, ma sempre su un solo core
//
const double K = 1.0;
const double D = 1.0;

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
        if(vector->data != NULL) free(vector->data);
        free(vector);
    }
}

double box_muller();
void istogrammi_stampa(struct Vector *counts, double minimo, double massimo, long long total_points, int binNum, const char *name);
void momenti stampa(double m1, double m2, double m3, double m4, const char* series_type, const char* filename);

int main(int argc, char *argv[]) {
    int maxTime_ensemble, intervallini, N_series;
    int seed_ensemble, n_bins;
    double x0;

    printf("--- Serial Ensemble Simulation ---\n");

    FILE *infile = fopen("input_ensemble_serial.dat", "r");

    if (infile == NULL) {
        perror("Err: non esiste 'input_ensemble_serial.dat'");
        return 1;
    }
    fscanf(infile, "%d %d %d %d %lf %d", &maxTime_ensemble, &intervallini, &N_series, &n_bins, &x0, &seed_ensemble);
    fclose(infile);

    long long int short_N_steps = (long long int)maxTime_ensemble * intervallini;
    if (short_N_steps == 0) short_N_steps = 1;
    double DT = 1.0 / intervallini;

    printf("\n--- Media d'insieme (Eseguita in Seriale) ---\n");
    printf("Generazione di %d serie totali da %.2e punti ciascuna.\n", N_series, (double)short_N_steps);

    if (seed_ensemble <= 0) srand((unsigned int)time(NULL));
    else srand((unsigned int)seed_ensemble);

    double global_sums[4] = {0.0, 0.0, 0.0, 0.0};
    struct Vector *global_hist_counts = newVector(n_bins);

    if (global_hist_counts == NULL) {
        fprintf(stderr, "ERRORE: Fallimento allocazione memoria per global_hist_counts.\n");
        return 1;
    }

    for(size_t i = 0; i < global_hist_counts->size; ++i) global_hist_counts->data[i] = 0;

    double est_min = -10.0, est_max = 10.0;

    for (int i = 0; i < N_series; i++) {
        struct Vector *current_series = newVector(short_N_steps);

        if (current_series == NULL) {
            fprintf(stderr, "ERRORE: Fallimento allocazione memoria per una serie.\n");
            return 1;
        }

        current_series->data[0] = x0;

        for (long long int j = 1; j < short_N_steps; j++) {
            double x = current_series->data[j-1];
            double h_x = (x >= 0) ? -K : K;
            double noise = sqrt(2.0 * D * DT) * box_muller();
            current_series->data[j] = x + h_x * DT + noise;
        }

        // Accumulate results directly into global variables
        double binSize = (est_max - est_min) / n_bins;
        for (long long int j = 0; j < short_N_steps; j++) {
            double x = current_series->data[j];
            global_sums[0] += x;
            global_sums[1] += x * x;
            global_sums[2] += x * x * x;
            global_sums[3] += x * x * x * x;

            int bin = (int)((x - est_min) / binSize);
            if (bin >= n_bins) bin = n_bins - 1;
            if (bin < 0) bin = 0;
            global_hist_counts->data[bin]++;
        }
        delVector(current_series);
    }

    printf("Inizio analisi ensemble...\n");
    long long global_total_points = (long long) N_series * short_N_steps;

    double m1 = global_sums[0] / global_total_points;
    double m2 = global_sums[1] / global_total_points;
    double m3 = global_sums[2] / global_total_points;
    double m4 = global_sums[3] / global_total_points;

    momenti stampa(m1, m2, m3, m4, "Ensemble (Serial)", "moments_ensemble_serial.dat");
    istogrammi_stampa(global_hist_counts, est_min, est_max, global_total_points, n_bins, "ensemble_globale_serial");

    delVector(global_hist_counts);
    printf("Fine analisi e scrittura file output.\n");

    return 0;
}


void momenti stampa(double m1, double m2, double m3, double m4, const char* series_type, const char* filename) {
    FILE *outfile = fopen(filename, "w");
    if (outfile == NULL) {
        perror("Errore durante l'apertura del file dei momenti");
        return;
    }
    double varianza = m2 - (m1 * m1);
    fprintf(outfile, "Momenti %s\n", series_type);
    fprintf(outfile, "Media <x>:\t\t%.6f\n", m1);
    fprintf(outfile, "Media Quadratica <x^2>:\t%.6f\n", m2);
    fprintf(outfile, "Momento Terzo <x^3>:\t%.6f\n", m3);
    fprintf(outfile, "Momento Quarto <x^4>:\t%.6f\n", m4);
    fprintf(outfile, "Varianza (<x^2>-<x>^2):\t%.6f\n\n", varianza);
    fclose(outfile);
}

void istogrammi_stampa(struct Vector *counts, double minimo, double massimo, long long total_points, int binNum, const char *name) {
    if (counts == NULL || total_points == 0) return;
    if (binNum == 0) {
        printf("ATTENZIONE: n_bins e' zero, impossibile generare l'istogramma.\n");
        return;
    }

    long long calculated_total_counts = 0;
    for(int i = 0; i < binNum; i++) {
        calculated_total_counts += (long long)round(counts->data[i]);
    }
    printf("Controllo Istogramma: Conteggi totali nei bin = %.2e\n", (double)calculated_total_counts);
    printf("Controllo Istogramma: Punti totali attesi    = %.2e\n", (double)total_points);
    if (calculated_total_counts != total_points) {
        printf("ATTENZIONE: I conteggi totali non corrispondono! E' possibile che alcuni punti siano caduti fuori dal range dell'istogramma.\n");
    }

    char filename[50];
    snprintf(filename, 50, "istogramma_%s.dat", name);
    double binSize = (massimo - minimo) / binNum;

    FILE *fp4 = fopen(filename, "w");
    if (fp4 == NULL) {
        perror("Errore durante l'apertura del file istogramma");
        return;
    }

    fprintf(fp4, "# Centro_Bin\tDensita_di_Probabilita\n");

    double total_normalized_area = 0.0;
    for(int i=0; i < binNum; i++) {
        double bin_center = minimo + (i + 0.5) * binSize;
        double normalized_height = (counts->data[i] / (double)total_points) / binSize;
        total_normalized_area += normalized_height * binSize;
        fprintf(fp4, "%f\t%f\n", bin_center, normalized_height);
    }

    printf("Controllo Istogramma 2: L'area totale normalizzata calcolata e' %f.\n", total_normalized_area);
    fclose(fp4);
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
