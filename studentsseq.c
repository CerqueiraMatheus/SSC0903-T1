#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <omp.h>

typedef struct _Stats {
    int min;
    int max;
    double med;
    double avg;
    double dev;
} Stats;

void print_stats(int min, int max, double med, double avg, double dev) {
    printf("menor: %d, maior: %d, mediana: %.2lf, média: %.2lf e DP: %.2lf\n", min, max, med, avg, dev);
}

int min(int *scores, int start, int end) {
    int min = scores[start];

    for (int i = start; i <= end; i++)
        if (scores[i] < min)
            min = scores[i];

    return min;
}

int max(int *scores, int start, int end) {
    int max = scores[start];

    for (int i = start; i <= end; i++)
        if (scores[i] > max)
            max = scores[i];

    return max;
}

double median(int *scores, int start, int end) {
    int buckets[101] = { 0 };
    for (int i = start; i <= end; i++)
        buckets[scores[i]]++;

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

double average(int *scores, int start, int end) {
    double avg = 0;

    int t = 1;
    for (int i = start; i <= end; i++) {
        avg += (scores[i] - avg) / t;
        
        t++;
    }

    return avg;
}

double std_dev(int *scores, int avg, int start, int end) {
    double res = 0;

    for (int i = start; i <= end; i++) {
        res += pow(scores[i] - avg, 2);
    }
    res = sqrt(res / (end - start + 1));

    return res;
}

typedef double (*medfn) (int *scores, int start, int end);

void get_stats(int *scores, Stats *stats, int curr_start, int n_items) {
    stats->min = min(scores, curr_start, curr_start + n_items - 1);
    stats->max = max(scores, curr_start, curr_start + n_items - 1);
    stats->med = median(scores, curr_start, curr_start + n_items - 1);
    stats->avg = average(scores, curr_start, curr_start + n_items - 1);
    stats->dev = std_dev(scores, stats->avg, curr_start, curr_start + n_items - 1);
}

int main(void) {
    int n_regions, n_cities, n_students, seed;
    
    if (scanf("%d %d %d %d", &n_regions, &n_cities, &n_students, &seed));

    srand(seed);
    
    int *scores = calloc(n_students * n_cities * n_regions, sizeof(int));
    for (int i = 0; i < n_regions * n_cities * n_students; i++) {
        scores[i] = rand() % 101;
    }

    Stats *stats = calloc(n_cities * n_regions + n_regions + 1, sizeof(Stats)); // Per city, per region and global stats
    
    int best_region = n_regions * n_cities;
    int best_city[2] = {0, 0};

    double start = omp_get_wtime();
    for (int i = 0; i < n_regions; i++) {
        for (int j = 0; j < n_cities; j++) {
            // Per city
            int curr_city = (i * n_cities * n_students) + (j * n_students);
            get_stats(scores, &(stats[i * n_cities + j]), curr_city, n_students);
            
            if (stats[i * n_cities + j].avg > stats[best_city[0] * n_cities + best_city[1]].avg) {
                best_city[0] = i;
                best_city[1] = j;
            }
        }
        // Per region
        int curr_region = (i * n_cities * n_students);
        get_stats(scores, &(stats[n_regions * n_cities + i]), curr_region, n_cities * n_students);

        if (stats[n_regions * n_cities + i].avg > stats[n_regions * n_cities + best_region].avg) {
            best_region = i;
        }
    }

    // Global
    get_stats(scores, &(stats[n_regions * n_cities + n_regions]), 0, n_regions * n_cities * n_students);
    double stop = omp_get_wtime();

    for (int i = 0; i < n_regions; i++) {
        for (int j = 0; j < n_cities; j++) {
            printf("Reg %d - Cid %d: ", i, j);
            print_stats(stats[i * n_cities + j].min, stats[i * n_cities + j].max, 
                        stats[i * n_cities + j].med, stats[i * n_cities + j].avg, 
                        stats[i * n_cities + j].dev);
        }
        printf("\n");
    }

    for (int i = 0; i < n_regions; i++) {
        if (stats[n_regions * n_cities + i].avg > stats[best_region].avg) {

        }
        printf("Reg %d: ", i);
        print_stats(stats[n_regions * n_cities + i].min, stats[n_regions * n_cities + i].max, 
                    stats[n_regions * n_cities + i].med, stats[n_regions * n_cities + i].avg, 
                    stats[n_regions * n_cities + i].dev);
    }
    printf("\n");

    printf("Brasil: ");
    print_stats(stats[n_regions * n_cities + n_regions].min, stats[n_regions * n_cities + n_regions].max, 
                stats[n_regions * n_cities + n_regions].med, stats[n_regions * n_cities + n_regions].avg, 
                stats[n_regions * n_cities + n_regions].dev);

    printf("\n");
    printf("Melhor regiao: Regiao %d\n", best_region);
    printf("Melhor cidade: Regiao %d, Cidade %d\n", best_city[0], best_city[1]);

    printf("\n");
    printf("Tempo de resposta sem considerar E/S, em segundos: %.4fs\n", (double)(stop-start));

    free(scores);
    free(stats);
}