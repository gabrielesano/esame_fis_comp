#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <mpi.h> // libreria per parallelizzazione

//
//Prova di concetto per la parallelizzazione
//Riguarda solo la generazioni delle serie di ensemble
//Ci sono diverse approssimazioni rispetto a quanto utilizzato per la risoluzione:
//per esempio i max e min per l'istogramma sono stimati, piuttosto che calcolati;
//possibile soluzione: Si potrebbe parallelizzare il controllo di massimo e minimo per ogni core,
//riunire i dati e dare massimi e minimi globali per poi successivamente effettuare un binning
//(Parallelo o no, da valutare l'overhead della comunicazione su questa operazione)
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

// prototipi
double box_muller();
void istogrammi_stampa(struct Vector *counts, double minimo, double massimo, long long total_points, int binNum, const char *name);
void momenti_stampa(double m1, double m2, double m3, double m4, const char* series_type, FILE *outfile);

int main(int argc, char *argv[]) {
    int rank, p;

    // Inizializzazione MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    // dichiarazioni
    int maxTime_ensemble, intervallini, N_series;
    int seed_ensemble, n_bins;
    double x0;

    // lettura input affidata al master
    if (rank == 0) {
        printf("Simulazione parallela su %d core\n", p);
        FILE *infile = fopen("input_ensemble.dat", "r");
        if (infile == NULL) {
            perror("Err: non esiste 'input_ensemble.dat'");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        fscanf(infile, "%d %d %d %d %lf %d", &maxTime_ensemble, &intervallini, &N_series, &n_bins, &x0, &seed_ensemble);
        fclose(infile);

        if (N_series % p != 0) {
            fprintf(stderr, "ERRORE: Il numero di serie (N_series = %d) deve essere un multiplo del numero dei core (%d).\n", N_series, p);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    //fuori da if rank=0 tutti i core stanno ricevendo queste istruzioni
    // Funzione che comunica i dati dal master a tutti i core, viene eseguita da tutti
    MPI_Bcast(&maxTime_ensemble, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&intervallini, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N_series, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n_bins, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&x0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&seed_ensemble, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // parametri per la simulazione
    long long int short_N_steps = (long long int)maxTime_ensemble * intervallini;
    if (short_N_steps == 0) short_N_steps = 1;
    double DT = 1.0 / intervallini;

    if (rank == 0) {
        printf("\nMedia d'insieme (Eseguita in Parallelo)\n");
        printf("Generazione di %d serie totali da %.2e punti ciascuna.\n", N_series, (double)short_N_steps);
    }

    int series_per_proc = N_series / p;
    if (seed_ensemble <= 0) srand((unsigned int)time(NULL) + rank);
    else srand((unsigned int)seed_ensemble + rank);


    //dichiaro un array locale che aggiorna i valori dei momenti su ogni core
    double local_sums[4] = {0.0, 0.0, 0.0, 0.0};
    struct Vector *local_hist_counts = newVector(n_bins);

    if (local_hist_counts == NULL) {
        fprintf(stderr, "ERRORE (Rank %d): errore allocazione memoria per istogramam.\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for(size_t i = 0; i < local_hist_counts->size; ++i) local_hist_counts->data[i] = 0;

    double est_min = -10.0, est_max = 10.0; //stimo massimo e minimo per l'istogramma su valori tipici

    //integro il processo con drift k su ogni core
    //i momenti vengono calcolati e aggiornati sul momento al calcolo di ogni serie
    //per non eccedere la memoria che è possibile allocare
    for (int i = 0; i < series_per_proc; i++) {
        struct Vector *current_series = newVector(short_N_steps);

        if (current_series == NULL) {
            fprintf(stderr, "ERRORE (Rank %d): Fallimento allocazione memoria per una serie.\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        current_series->data[0] = x0;

        for (long long int j = 1; j < short_N_steps; j++) {
            double x = current_series->data[j-1];
            double h_x = (x >= 0) ? -K : K;
            double noise = sqrt(2.0 * D * DT) * box_muller();
            current_series->data[j] = x + h_x * DT + noise;
        }

        double binSize = (est_max - est_min) / n_bins;
        for (long long int j = 0; j < short_N_steps; j++) {
            double x = current_series->data[j];
            local_sums[0] += x;
            local_sums[1] += x * x;
            local_sums[2] += x * x * x;
            local_sums[3] += x * x * x * x;

            int bin = (int)((x - est_min) / binSize);
            if (bin >= n_bins) bin = n_bins - 1;
            if (bin < 0) bin = 0;
            local_hist_counts->data[bin]++;
        }
        delVector(current_series);
    }

    // Riduzione dei valori calcolati su ogni core
    // qui ci sono le operazioni effettuate dal singolo core
    double global_sums[4];
    struct Vector *global_hist_counts = NULL;
    if (rank == 0) {
        global_hist_counts = newVector(n_bins);

        if (global_hist_counts == NULL) {
            fprintf(stderr, "ERRORE (Rank %d): errore allocazione memoria per global_hist_counts.\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    MPI_Reduce(local_sums, global_sums, 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    //ricezione conteggi bin istogramma
    double *recv_buffer = NULL;
    if (rank == 0) {
        recv_buffer = global_hist_counts->data;
    }
    MPI_Reduce(local_hist_counts->data, recv_buffer, n_bins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    delVector(local_hist_counts);

    // Il core 0 va a mediare tutti i dati e stampare l'istogramma totale
    if (rank == 0) {
        printf("Inizio analisi ensemble...\n");
        long long global_total_points = (long long) N_series * short_N_steps;

        double m1 = global_sums[0] / global_total_points;
        double m2 = global_sums[1] / global_total_points;
        double m3 = global_sums[2] / global_total_points;
        double m4 = global_sums[3] / global_total_points;

        FILE *moments_ensemble_file = fopen("moments_ensemble.dat", "w");
        momenti_stampa(m1, m2, m3, m4, "Ensemble", moments_ensemble_file);
        fclose(moments_ensemble_file);

        istogrammi_stampa(global_hist_counts, est_min, est_max, global_total_points, n_bins, "ensemble_globale");

        delVector(global_hist_counts);
        printf("Fine analisi e scrittura file output.\n");
    }

    // --- Finalize MPI ---
    MPI_Finalize();
    return 0;
}

// --- Helper Functions ---

void momenti_stampa(double m1, double m2, double m3, double m4, const char* series_type, FILE *outfile) {
    if (outfile == NULL) return;
    double varianza = m2 - (m1 * m1);
    fprintf(outfile, "Momenti %s\n", series_type);
    fprintf(outfile, "Media <x>:\t\t%.6f\n", m1);
    fprintf(outfile, "Media Quadratica <x^2>:\t%.6f\n", m2);
    fprintf(outfile, "Momento Terzo <x^3>:\t%.6f\n", m3);
    fprintf(outfile, "Momento Quarto <x^4>:\t%.6f\n", m4);
    fprintf(outfile, "Varianza (<x^2>-<x>^2):\t%.6f\n\n", varianza);
}

void istogrammi_stampa(struct Vector *counts, double minimo, double massimo, long long total_points, int binNum, const char *name) {
    if (counts == NULL || total_points == 0) return;
    // Add a safety check to prevent division by zero if n_bins is 0
    if (binNum == 0) {
        printf("ATTENZIONE: n_bins e' zero, impossibile generare l'istogramma.\n");
        return;
    }

    // controllo punti totali nei bin
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

    // controllo normalizzazione area
    double total_normalized_area = 0.0;

    for(int i=0; i < binNum; i++) {
        double bin_center = minimo + (i + 0.5) * binSize;
        double normalized_height = (counts->data[i] / (double)total_points) / binSize;

        // Accumulate the area of each bin (height * width)
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
