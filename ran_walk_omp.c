#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#define ACCURACY 1000
#include <assert.h>
#include <sys/time.h>


int main(int argc, char ** argv) {
    int a, b, x, N, P, b_count = 0;
    double p, chance_b, start_time, time, step = 0;

    omp_lock_t lock;

    if (argc != 7)
    {
        printf("Wrong number of args\n");
        return 0;
    }
    else
    {
        a = atoi(argv[1]);
        b = atoi(argv[2]);
        x = atoi(argv[3]);
        N = atoi(argv[4]);
        p = atof(argv[5]);
        P = atoi(argv[6]);
    }

    int * rands = (int*)malloc(sizeof(int) * P);
    for (size_t i = 0; i < P; i++)
    {
        rands[i] = rand();
    }

    start_time = omp_get_wtime();
    omp_init_lock(&lock);
    omp_set_num_threads(P);

    #pragma omp parallel for

    for(unsigned int i = 0; i < N; i++)
    {
        int curr_pos = x;
        unsigned int next_step;
        unsigned int steps_done = 0;
        unsigned int seed = (unsigned int) clock() * rands[omp_get_thread_num()];

     //захватываем lock только N раз, когда каждая частица завершила свое блуждание

        while (curr_pos != b && curr_pos != a)
        {
            next_step = (unsigned) rand_r(&seed) % ACCURACY < ACCURACY * p ? 0 : 1;
            if (next_step == 0) {
                curr_pos++;
            } else {
                curr_pos--;
            }
            steps_done++;
        }

        if (curr_pos == b)
        {
            omp_set_lock(&lock);
            chance_b++;
            step += steps_done;
            omp_unset_lock(&lock);
        }

    }

    time = omp_get_wtime() - start_time;

    chance_b = (double) b_count / (double) N;
    step /= (double) N;

    FILE * out = fopen("stats.txt", "w");
    if (out == NULL)
    {
        printf("Can't open the file\n");
        exit(1);
    }
    fprintf(out, "%.2f %.1lf %lfs %d %d %d %d %.2f %d\n", chance_b, step, time, a, b, x, N, p, P);
    fclose(out);

    printf("total time: %lf\n", time);

    free(rands);

    return 0;
}
