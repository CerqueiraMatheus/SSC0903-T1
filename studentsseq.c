#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include<limits.h>

typedef struct _Stats {
    int min;
    int max;
    double med;
    double avg;
    double dev;
} Stats;

typedef struct _HeapTuple {
    int data;
    int list_idx;
} HeapTuple;

typedef struct _Heap {
    HeapTuple *array;
    int size;
} Heap;

int compare_tuple(HeapTuple *a, HeapTuple *b);

void merge_sorted_lists(int *scores, int start, int n_lists, int n_elems);
void merge_sort(int *scores, int left, int mid, int right);
void partition(int *scores, int left, int right);

Heap *heap_init(int heapsize);
void heap_free(Heap *heap);
void heap_insert(Heap *h, HeapTuple element);
void heap_pop(Heap* h, HeapTuple *ht_res);


void print_stats(int min, int max, double med, double avg, double dev) {
    printf("menor: %d, maior: %d, mediana: %.2lf, média: %.2lf e DP: %.2lf\n", min, max, med, avg, dev);
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

void get_stats(int *scores, Stats *stats, medfn medfunc, int curr_start, int n_items) {
    partition(scores, curr_start, curr_start + n_items - 1);

    stats->min = scores[curr_start];
    stats->max = scores[curr_start + n_items - 1];
    stats->med = medfunc(scores, curr_start, curr_start + n_items - 1);
    stats->avg = average(scores, curr_start, curr_start + n_items - 1);
    stats->dev = std_dev(scores, stats->avg, curr_start, curr_start + n_items - 1);
}

// int main(void) {
//     int n_regions, n_cities, n_students, seed;
    
//     scanf("%d %d %d %d", &n_regions, &n_cities, &n_students, &seed);

//     srand(seed);
    
//     int scores[n_regions * n_cities * n_students];
//     for (int i = 0; i < n_regions * n_cities * n_students; i++) {
//         scores[i] = rand() % 101;
//     }

//     medfn medfunc_city = n_students % 2 ? median_odd : median_even;
//     medfn medfunc_region = n_cities % 2 ? medfunc_city : median_even;
//     medfn medfunc_global = n_regions % 2 ? medfunc_region : median_even;

//     Stats stats[1 + n_regions + n_cities * n_regions ]; // Per city, per region and global stats
    
//     int best_region = n_regions * n_cities;
//     int best_city[2] = {0, 0};

//     clock_t start=clock();
//     for (int i = 0; i < n_regions; i++) {
//         for (int j = 0; j < n_cities; j++) {
//             // Per city
//             int curr_city = (i * n_cities * n_students) + (j * n_students);
//             get_stats(scores, &(stats[i * n_cities + j]), medfunc_city, curr_city, n_students);
            
//             if (stats[i * n_cities + j].avg > stats[best_city[0] * n_cities + best_city[1]].avg) {
//                 best_city[0] = i;
//                 best_city[1] = j;
//             }
//         }
//         // Per region
//         int curr_region = (i * n_cities * n_students);
//         get_stats(scores, &(stats[n_regions * n_cities + i]), medfunc_region, curr_region, n_cities * n_students);

//         if (stats[n_regions * n_cities + i].avg > stats[n_regions * n_cities + best_region].avg) {
//             best_region = i;
//         }
//     }

//     // Global
//     get_stats(scores, &(stats[n_regions * n_cities + n_regions]), medfunc_global, 0, n_regions * n_cities * n_students);
//     clock_t stop = clock();

//     for (int i = 0; i < n_regions; i++) {
//         for (int j = 0; j < n_cities; j++) {
//             printf("Reg %d - Cid %d: ", i, j);
//             print_stats(stats[i * n_cities + j].min, stats[i * n_cities + j].max, 
//                         stats[i * n_cities + j].med, stats[i * n_cities + j].avg, 
//                         stats[i * n_cities + j].dev);
//         }
//         printf("\n");
//     }

//     for (int i = 0; i < n_regions; i++) {
//         if (stats[n_regions * n_cities + i].avg > stats[best_region].avg) {

//         }
//         printf("Reg %d: ", i);
//         print_stats(stats[n_regions * n_cities + i].min, stats[n_regions * n_cities + i].max, 
//                     stats[n_regions * n_cities + i].med, stats[n_regions * n_cities + i].avg, 
//                     stats[n_regions * n_cities + i].dev);
//     }
//     printf("\n");

//     printf("Brasil: ");
//     print_stats(stats[n_regions * n_cities + n_regions].min, stats[n_regions * n_cities + n_regions].max, 
//                 stats[n_regions * n_cities + n_regions].med, stats[n_regions * n_cities + n_regions].avg, 
//                 stats[n_regions * n_cities + n_regions].dev);

//     printf("\n");
//     printf("Melhor regiao: Regiao %d\n", best_region);
//     printf("Melhor cidade: Regiao %d, Cidade %d\n", best_city[0], best_city[1]);

//     printf("\n");
//     printf("Tempo de resposta sem considerar E/S, em segundos: %.4fs\n", (double)(stop-start)/CLOCKS_PER_SEC);
// }

int main(void) {
    int scores[12] = {1, 3, 5, 7, 2, 8, 10, 11, 4, 6, 9, 12};
    merge_sorted_lists(scores, 0, 3, 4);
    
    printf("\n\n");
    for (int i = 0; i < 12; i++)
        printf("%d ", scores[i]);
    printf("\n\n");
}

void merge_sorted_lists(int *scores, int start, int n_lists, int n_elems) {
    Heap *heap = heap_init(n_lists);
    int *curr_idx_lists = calloc(n_lists, sizeof(int));

    // [0...n_elems - 1]
    HeapTuple aux;
    for (int i = 0; i < n_lists; i++) {
        aux.data = scores[start + i * n_elems];
        aux.list_idx = i;
        heap_insert(heap, aux);
    }

    int sorted_scores[n_lists * n_elems];
    int curr_pos = 0;
    while (heap->size > 0) {
        heap_pop(heap, &aux);
        sorted_scores[curr_pos++] = aux.data;

        curr_idx_lists[aux.list_idx]++;
        if (curr_idx_lists[aux.list_idx] < n_elems - 1) {
            aux.data = curr_idx_lists[aux.list_idx];
            heap_insert(heap, aux);
        }
    }

    heap_free(heap);
    free(curr_idx_lists);

    for (int i = 0; i < n_lists * n_elems; i++) {
        scores[start + i] = sorted_scores[i];
    }
}

void merge_sort(int scores[], int left, int mid, int right) {
    int i, j, k;
    int n1 = mid - left + 1;
    int n2 = right - mid;
 
    int L[n1], R[n2];
 
    for (i = 0; i < n1; i++)
        L[i] = scores[left + i];
    for (j = 0; j < n2; j++)
        R[j] = scores[mid + 1 + j];
 
    i = 0;
    j = 0;
    k = left;
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
 
    while (i < n1) {
        scores[k] = L[i];
        i++;
        k++;
    }
 
    while (j < n2) {
        scores[k] = R[j];
        j++;
        k++;
    }
}
 
void partition(int scores[], int left, int right) {
    if (left < right) {
        int m = left + (right - left) / 2;
 
        partition(scores, left, m);
        partition(scores, m + 1, right);
 
        merge_sort(scores, left, m, right);
    }
}

/*Initialize Heap*/
Heap *heap_init(int heapsize) {
    Heap *heap = calloc(1, sizeof(Heap));

    heap->array = calloc(heapsize + 1, sizeof(int));
    heap->array[0].data = -INT_MAX;
    heap->array[0].list_idx = -INT_MAX;

    heap->size = 0;

    return heap;
}

void heap_free(Heap *heap) {
    free(heap->array);
    free(heap);
}
 
/*Insert an element into the heap */
void heap_insert(Heap *h, HeapTuple ht) {
    h->size++;
    // ! FALTA CHECAGEM HEAP TAMANHO CORRETO

    h->array[h->size].data = ht.data ; /*Insert in the last place*/
    h->array[h->size].list_idx = ht.list_idx; 
    /*Adjust its position*/
    int now = h->size;
    while (h->array[now / 2].data > ht.data) {
        h->array[now].data = h->array[now / 2].data;
        h->array[now].list_idx = h->array[now / 2].list_idx;
        now /= 2;
    }
    h->array[now].data = ht.data;
    h->array[now].list_idx = ht.list_idx;
}
 
void heap_pop(Heap* h, HeapTuple *ht_res) {
    /* heap[1] is the minimum element. So we remove heap[1]. Size of the heap is decreased.
     Now heap[1] has to be filled. We put the last element in its place and see if it fits.
     If it does not fit, take minimum element among both its children and replaces parent with it.
     Again See if the last element fits in that place.*/
    HeapTuple min_element, last_element;
    int child, now;

    min_element = h->array[1];
    last_element = h->array[h->size--];
    /* now refers to the index at which we are now */
    for (now = 1; now * 2 <= h->size; now = child) {
        /* child is the index of the element which is minimum among both the children */
        /* Indexes of children are i*2 and i*2 + 1*/
        child = now * 2;
        /*child!=heapSize beacuse heap[heapSize+1] does not exist, which means it has only one
         child */
        if (child != h->size && h->array[child + 1].data < h->array[child].data)
            child++;
        
        /* To check if the last element fits ot not it suffices to check if the last element
         is less than the minimum element among both the children*/
        if (last_element.data > h->array[child].data) {
            h->array[now].data = h->array[child].data;
            h->array[now].list_idx = h->array[child].list_idx;
        }
        else 
            break;
        
    }
    h->array[now].data = last_element.data;
    h->array[now].list_idx = last_element.list_idx;

    ht_res->data = min_element.data;
    ht_res->list_idx = min_element.list_idx;
}

int compare_tuple(HeapTuple *a, HeapTuple *b) {
    return a->data - b->data;
}