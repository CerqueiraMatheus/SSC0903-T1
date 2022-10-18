#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void merge(int *scores, int left, int mid, int right);
void mergeSort(int *scores, int left, int right);

typedef struct _Stats {
    int min;
    int max;
    double med;
    double avg;
    double dev;
} Stats;

typedef struct _TopStats {
    int best_region;
    int best_city[2];
    double best_region_avg;
    double best_city_avg;
} TopStats;


void print_stats(int min, int max, double med, double avg, double dev) {
    printf("menor: %d, maior: %d, mediana: %.2lf, média: %.2lf e DP: %.2lf\n", min, max, med, avg, dev);
}

int cmpint (const void *a, const void *b) {
   const int *A = a, *B = b;
   return *A - *B;
}

double median_even(int *scores, int start, int end) {
    int med = (end - start + 1) / 2;
    return (double) (scores[start + med] + scores[start + med - 1]) / 2;
}

double median_odd(int *scores, int start, int end) {
    return scores[start + (int) (end - start + 1) / 2];
}

double average(int *scores, int start, int end) {
    int summ = 0;
    for (int i = start; i <= end; i++)
        summ += scores[i];
    
    return (double) summ / (end - start + 1);
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

int get_stats(int *scores, Stats *stats, medfn medfunc, int curr_start, int n_items) {
    mergeSort(scores, curr_start, curr_start + n_items - 1);

    stats->min = scores[curr_start];
    stats->max = scores[curr_start + n_items - 1];
    stats->med = medfunc(scores, curr_start, curr_start + n_items - 1);
    stats->avg = average(scores, curr_start, curr_start + n_items - 1);
    stats->dev = std_dev(scores, stats->avg, curr_start, curr_start + n_items - 1);
}

int main(void) {
    int n_regions, n_cities, n_students, seed;
    
    scanf("%d %d %d %d", &n_regions, &n_cities, &n_students, &seed);

    int scores[n_regions * n_cities * n_students];

    srand(seed);
    for (int i = 0; i < n_regions * n_cities * n_students; i++) {
        scores[i] = rand() % 100;
    }

    medfn medfunc_city = n_students % 2 ? median_odd : median_even;
    medfn medfunc_region = n_cities % 2 ? medfunc_city : median_even;
    medfn medfunc_global = n_regions % 2 ? medfunc_region : median_even;

    Stats stats[n_regions + 2]; // Local, global and per region
    TopStats top_stats;

    for (int i = 0; i < n_regions; i++) {
        
        for (int j = 0; j < n_cities; j++) {
            int curr_city = (i * n_cities * n_students) + (j * n_students);
            get_stats(scores, &(stats[0]), medfunc_city, curr_city, n_students);

            printf("Reg %d - Cid %d: ", i, j);
            print_stats(stats[0].min, stats[0].max, stats[0].med, stats[0].avg, stats[0].dev);
            
            if (top_stats.best_city_avg < stats[0].avg) {
                top_stats.best_city[0] = i;
                top_stats.best_city[1] = j;
            }
        }
        printf("\n");

        // Per region
        int curr_region = (i * n_cities * n_students);
        get_stats(scores, &(stats[i + 2]), medfunc_region, curr_region, n_cities * n_students);

        if (top_stats.best_region_avg < stats[i + 2].avg) {
            top_stats.best_region = i;
        }
    }

    // Global
    get_stats(scores, &(stats[1]), medfunc_global, 0, n_regions * n_cities * n_students);

    for (int i = 0; i < n_regions; i++) {
        printf("Reg %d: ", i);
        print_stats(stats[i + 2].min, stats[i + 2].max, stats[i + 2].med, stats[i + 2].avg, stats[i + 2].dev);
    }
    printf("\n");

    printf("Brasil: ");
    print_stats(stats[1].min, stats[1].max, stats[1].med, stats[1].avg, stats[1].dev);

    printf("\n");
    printf("Melhor regiao: Regiao %d\n", top_stats.best_region);
    printf("Melhor cidade: Regiao %d, Cidade %d\n", top_stats.best_city[0], top_stats.best_city[1]);

    printf("Tempo de resposta: xablau");
}

void merge(int scores[], int left, int mid, int right) {
    int i, j, k;
    int n1 = mid - left + 1;
    int n2 = right - mid;
 
    /* create temp arrays */
    int L[n1], R[n2];
 
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++)
        L[i] = scores[left + i];
    for (j = 0; j < n2; j++)
        R[j] = scores[mid + 1 + j];
 
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = left; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            scores[k] = L[i];
            i++;
        }
        else {
            scores[k] = R[j];
            j++;
        }
        k++;
    }
 
    /* Copy the remaining elements of L[], if there
    are any */
    while (i < n1) {
        scores[k] = L[i];
        i++;
        k++;
    }
 
    /* Copy the remaining elements of R[], if there
    are any */
    while (j < n2) {
        scores[k] = R[j];
        j++;
        k++;
    }
}
 
/* l is for left index and r is right index of the
sub-array of arr to be sorted */
void mergeSort(int scores[], int left, int right) {
    if (left < right) {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        int m = left + (right - left) / 2;
 
        // Sort first and second halves
        mergeSort(scores, left, m);
        mergeSort(scores, m + 1, right);
 
        merge(scores, left, m, right);
    }
}