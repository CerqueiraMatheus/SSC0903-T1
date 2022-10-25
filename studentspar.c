#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <omp.h>

#define N_SCORES 101

typedef struct _Stats {
    int n;
    int min;
    int max;
    double med;
    double avg;
    double var;
} Stats;

#pragma omp declare reduction(reduce_stats: Stats: \
    omp_out.min = omp_out.min <= omp_in.min ? omp_out.min : omp_in.min, \
    omp_out.max = omp_out.max >= omp_in.max ? omp_out.max : omp_in.max, \
    omp_out.var = ((omp_out.var * omp_out.n + omp_in.var * omp_in.n) / (omp_out.n + omp_in.n)) + ((omp_out.n * omp_in.n * ((omp_out.avg - omp_in.avg)*(omp_out.avg - omp_in.avg))) / ((omp_out.n + omp_in.n)*(omp_out.n + omp_in.n))), \
    omp_out.avg = (omp_out.avg + omp_in.avg) / 2, \
    omp_out.n += omp_in.n) \
    initializer( omp_priv = { 0, INT_MAX, INT_MIN, 0, 0, 0 })

void print_stats(Stats stats) {
    printf("menor: %d, maior: %d, mediana: %.2lf, média: %.2lf e DP: %.2lf\n", stats.min, stats.max, stats.med, stats.avg, sqrt(stats.var));
}

int *counting_sort(int *scores, int* buckets, int start, int end) {
    for (int i = start; i <= end; i++)
        buckets[scores[i]]++;

    return buckets;
}

int get_min(int *buckets) {
    int i = 0;
    for (i = 0; i < N_SCORES; i++)
        if (buckets[i] > 0)
            break;
    
    return i;
}

int get_max(int *buckets) {
    int i = 0;
    for (i = N_SCORES - 1; i >= 0; i--)
        if (buckets[i] > 0)
            break;

    return i;
}

double get_median(int *buckets, int n_elems) {
    int start = 0;
    int end = n_elems - 1;
    int left = -1, right = -1;

    int curr_idx = -1;
    for (int i = 0; i <= 100; i++) {
        curr_idx += buckets[i];

        if (left == -1 && curr_idx >= (end - start) / 2) {
            left = i;
        }
        if (curr_idx >= (end - start + 1) / 2) {
            right = i;
            break;
        }
    }

    return (double) (left + right) / 2;
}

double get_avg(int *buckets, int n_elems) {
    int summ = 0;
    for (int i = 0; i < N_SCORES; i++)
        summ += i * buckets[i];
    
    return (double) summ / n_elems;
}

double get_var(int *buckets, int avg, int n_elems) {
    double res = 0;

    for (int i = 0; i < N_SCORES; i++) {
        res += buckets[i] * pow(i - avg, 2);
    }
    return res / n_elems;
}

typedef double (*medfn) (int *buckets, int n_elems);

int main(void) {
    int n_regions, n_cities, n_students, seed;
    
    scanf("%d %d %d %d", &n_regions, &n_cities, &n_students, &seed);

    srand(seed);
    
    int scores[n_regions * n_cities * n_students];
    for (int i = 0; i < n_regions * n_cities * n_students; i++) {
        scores[i] = rand() % 101;

        printf("%d - ", scores[i]);

        if ((i + 1) % n_students == 0)
            printf("\n");
    }

    Stats *stats = calloc(n_cities * n_regions + n_regions + 1, sizeof(Stats));

    int best_region = n_regions * n_cities;
    int best_city[2] = {0, 0};

    clock_t start=clock();

    int buckets[N_SCORES] = { 0 };
    Stats curr_stats;


    #pragma omp parallel for \
            reduction(reduce_stats: curr_stats) reduction(+: buckets)
    for (int i = 0; i < n_regions; i++) {

        #pragma omp parallel for \
            reduction(reduce_stats: curr_stats) reduction(+: buckets)
        for (int j = 0; j < n_cities; j++) {
            int curr_city = (i * n_cities * n_students) + (j * n_students);

            counting_sort(scores, buckets, curr_city, curr_city + n_students - 1);

            double avg = get_avg(buckets, n_students);
            curr_stats = (Stats){ 
                .n = n_students, 
                .min = get_min(buckets), 
                .max = get_max(buckets), 
                .med = get_median(buckets, n_students), 
                .avg = avg,
                .var = get_var(buckets, avg, n_students) 
            };
            printf("Reg %d City %i - ", i, j);
            print_stats(curr_stats);

            stats[i * n_cities + j] = curr_stats;
        }

        curr_stats.med = get_median(buckets, n_cities * n_students);
        stats[n_regions * n_cities + i] = curr_stats;
    } // End of parallel section

    curr_stats.med = get_median(buckets, n_regions * n_cities * n_students);
    stats[n_regions * n_cities + n_regions] = curr_stats;

    clock_t stop = clock();

    printf("\n");
    printf("\n");


    for (int i = 0; i < n_regions; i++) {
        for (int j = 0; j < n_cities; j++) {
            printf("Reg %d - Cid %d: ", i, j);
            print_stats(stats[i * n_cities + j]);
        }
        printf("\n");
    }

    for (int i = 0; i < n_regions; i++) {
        if (stats[n_regions * n_cities + i].avg > stats[best_region].avg) {

        }
        printf("Reg %d: ", i);
        print_stats(stats[n_regions * n_cities + i]);
    }
    printf("\n");

    printf("Brasil: ");
    print_stats(stats[n_regions * n_cities + n_regions]);

    printf("\n");
    printf("Melhor regiao: Regiao %d\n", best_region);
    printf("Melhor cidade: Regiao %d, Cidade %d\n", best_city[0], best_city[1]);

    printf("\n");
    printf("Tempo de resposta sem considerar E/S, em segundos: %.4fs\n", (double)(stop-start)/CLOCKS_PER_SEC);

    free(stats);
}