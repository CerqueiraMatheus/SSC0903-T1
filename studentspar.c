#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <omp.h>

#define N_SCORES 101

typedef struct _Stats {
    double n;
    int min;
    int max;
    double med;
    double avg;
    double var;
} Stats;

typedef struct Avg {
    int pos;
    double avg;
} Avg;

#pragma omp declare reduction(best: Avg: \
    omp_out.avg = (omp_in.avg >= omp_out.avg ? omp_in.avg : omp_out.avg), \
    omp_out.pos = (omp_in.avg >= omp_out.avg ? omp_in.pos : omp_out.pos)) \
    initializer( omp_priv = { 0, 0 })

#pragma omp declare reduction(reduce_stats: Stats: \
    omp_out.min = (omp_out.min <= omp_in.min ? omp_out.min : omp_in.min), \
    omp_out.max = (omp_out.max >= omp_in.max ? omp_out.max : omp_in.max), \
    omp_out.var = omp_out.n == 0 ? omp_in.var : ((omp_out.var * omp_out.n + omp_in.var * omp_in.n) / (omp_out.n + omp_in.n)) + ((omp_out.n * omp_in.n * ((omp_out.avg - omp_in.avg)*(omp_out.avg - omp_in.avg))) / ((omp_out.n + omp_in.n)*(omp_out.n + omp_in.n))), \
    omp_out.avg = omp_out.n == 0 ? omp_in.avg : (omp_out.n * omp_out.avg + omp_in.n * omp_in.avg) / (omp_out.n + omp_in.n), \
    omp_out.n = omp_out.n + omp_in.n) \
    initializer( omp_priv = { 0, INT_MAX, INT_MIN, 0, 0, 0 })

void reduce_stats(Stats* regional_stats, Stats* local_stats) {
    regional_stats->min = (regional_stats->min <= local_stats->min ? regional_stats->min : local_stats->min);
    regional_stats->max = (regional_stats->max >= local_stats->max ? regional_stats->max : local_stats->max);
    regional_stats->var = regional_stats->n == 0 ? local_stats->var : ((regional_stats->var * regional_stats->n + local_stats->var * local_stats->n) / (regional_stats->n + local_stats->n)) + ((regional_stats->n * local_stats->n * ((regional_stats->avg - local_stats->avg)*(regional_stats->avg - local_stats->avg))) / ((regional_stats->n + local_stats->n)*(regional_stats->n + local_stats->n)));
    regional_stats->avg = regional_stats->n == 0 ? local_stats->avg : (regional_stats->n * regional_stats->avg + local_stats->n * local_stats->avg) / (regional_stats->n + local_stats->n);
    regional_stats->n = regional_stats->n + local_stats->n;
}

void print_stats(Stats stats) {
    printf("menor: %d, maior: %d, mediana: %.2lf, m??dia: %.2lf e DP: %.2lf\n", stats.min, stats.max, stats.med, stats.avg, sqrt(stats.var));
}

void counting_sort(int *scores, int* buckets, int* local_buckets, int start, int end) {
    for (int i = start; i <= end; i++) {
        buckets[scores[i]]++;
        local_buckets[scores[i]]++;
    }
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

void print_buckets(int *buckets) {
    for (int i = 0; i < N_SCORES; i++)
        printf("b[%d]=%d - ", i, buckets[i]);

    printf("\n");
}

typedef double (*medfn) (int *buckets, int n_elems);

int main(void) {
    int n_regions, n_cities, n_students, seed;
    
    if (scanf("%d %d %d %d", &n_regions, &n_cities, &n_students, &seed));

    srand(seed);
    
    int *scores = calloc(n_students * n_cities * n_regions, sizeof(int));
    for (int i = 0; i < n_regions * n_cities * n_students; i++) {
        scores[i] = rand() % 101;
    }

    Stats *stats = calloc(n_cities * n_regions + n_regions + 1, sizeof(Stats));
    Avg best_city = { 0, 0 }, best_region = { 0, 0 };

    int global_buckets[N_SCORES] = { 0 };
    Stats global_stats = { 0, INT_MAX, INT_MIN, 0, 0, 0 };

    double start= omp_get_wtime();
    omp_set_nested(2);

    #pragma omp parallel for \
            default(none) \
            shared(n_regions, n_cities, n_students, scores, stats) \
            reduction(reduce_stats: global_stats) \
            reduction(+: global_buckets) \
            reduction(best: best_region) \
            reduction(best: best_city) \
            num_threads(4)
    for (int i = 0; i < n_regions; i++) {
        int regional_buckets[N_SCORES] = { 0 };
        Stats regional_stats = { 0, INT_MAX, INT_MIN, 0, 0, 0 };
        
        #pragma omp parallel for \
            default(none) \
            shared(n_regions, n_cities, n_students, scores, stats) \
            firstprivate(i) \
            reduction(reduce_stats: regional_stats) \
            reduction(+: regional_buckets) \
            reduction(best: best_city) \
            num_threads(2)
        for (int j = 0; j < n_cities; j++) {
            int local_buckets[N_SCORES] = { 0 };

            int curr_city = (i * n_cities * n_students) + (j * n_students);

            counting_sort(scores, regional_buckets, local_buckets, curr_city, curr_city + n_students - 1);

            double avg = get_avg(local_buckets, n_students);
            Stats local_stats = (Stats){ 
                .n = n_students, 
                .min = get_min(local_buckets),
                .max = get_max(local_buckets), 
                .med = get_median(local_buckets, n_students), 
                .avg = avg,
                .var = get_var(local_buckets, avg, n_students) 
            };

            reduce_stats(&regional_stats, &local_stats);
            
            if (local_stats.avg > best_city.avg) {
                best_city.avg = local_stats.avg;
                best_city.pos = i * n_cities + j;
            }

            stats[i * n_cities + j] = local_stats;
        }

        reduce_stats(&global_stats, &regional_stats);

        regional_stats.med = get_median(regional_buckets, n_cities * n_students);
        stats[n_regions * n_cities + i] = regional_stats;
        
        for (int i = 0; i < N_SCORES; i++)
            global_buckets[i] += regional_buckets[i];

        if (regional_stats.avg > best_region.avg) {
            best_region.avg = regional_stats.avg;
            best_region.pos = i;
        }
    } // End of parallel section

    global_stats.med = get_median(global_buckets, n_regions * n_cities * n_students);
    stats[n_regions * n_cities + n_regions] = global_stats;

    double stop = omp_get_wtime();

    printf("\n\n");

    for (int i = 0; i < n_regions; i++) {
        for (int j = 0; j < n_cities; j++) {
            printf("Reg %d - Cid %d: ", i, j);
            print_stats(stats[i * n_cities + j]);
        }
        printf("\n");
    }

    for (int i = 0; i < n_regions; i++) {
        printf("Reg %d: ", i);
        print_stats(stats[n_regions * n_cities + i]);
    }
    printf("\n");

    printf("Brasil: ");
    print_stats(stats[n_regions * n_cities + n_regions]);

    printf("\n");
    printf("Melhor regiao: Regiao %d\n", best_region.pos);
    printf("Melhor cidade: Regiao %d, Cidade %d\n", best_city.pos / n_cities, best_city.pos % n_cities);

    printf("\n");
    printf("Tempo de resposta sem considerar E/S, em segundos: %.4fs\n", (double)(stop-start));

    free(stats);
    free(scores);
}

