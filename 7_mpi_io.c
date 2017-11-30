#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

struct particle {
    int x;
    int y;
    int r;
};

struct particle creater_particle(int x, int y, int r) {
    struct particle arg;
    arg.x = x;
    arg.y = y;
    arg.r = r;
    return arg;
}

int main(int argc, char * argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int l = atoi(argv[1]);
    int a = atoi(argv[2]);
    int b = atoi(argv[3]);
    int N = atoi(argv[4]);

    double start = MPI_Wtime();
    int new_size = ((rank % a) * l + (rank / a) * a * l * l) * size;
    struct particle* particles = (struct particle*) malloc(sizeof(struct particle) * N);

    int seed;
    int* seeds = (int*)malloc(sizeof(int)* size);
    if (rank == 0) {
        srand(time(NULL));
        for (int i = 0; i < size; ++i) {
            seeds[i] = rand();
        }
    }

    MPI_Scatter(seeds, 1, MPI_UNSIGNED, &seed, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    free(seeds);

    for(int i = 0; i < N; ++i) {
        particles[i] = creater_particle(rand_r((unsigned int*)&seed) % l, rand_r((unsigned int*)&seed) % l,
                                  rand_r((unsigned int*)&seed) % (a * b));
    }

    int* array = (int*)malloc(sizeof(int) * l * l * size);
    for (int i = 0; i < l * l * size; i++) {
        array[i] = 0;
    }

    MPI_File file_bin;
    MPI_File_delete("data.bin", MPI_INFO_NULL);
    MPI_File_open(MPI_COMM_WORLD, "data.bin", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file_bin);

    for(int i = 0; i < N; ++i) {
        array[particles[i].y * l * size + particles[i].x * size + particles[i].r] += 1;
    }

    MPI_Aint intex;
    MPI_Aint e;
    MPI_Type_get_extent(MPI_INT, &e, &intex);
    MPI_Datatype view;
    MPI_Type_vector(l, l * size, l * a * size, MPI_INT, &view);
    MPI_Type_commit(&view);
    MPI_File_set_view(file_bin, new_size * sizeof(int), MPI_INT, view, "native", MPI_INFO_NULL);
    MPI_File_write(file_bin, array, l * l * size, MPI_INT, MPI_STATUS_IGNORE);
    MPI_Type_free(&view);
    MPI_File_close(&file_bin);

    free(array);

    double end = MPI_Wtime();

    if (rank == 0) {
        FILE *f_stats;
        f_stats = fopen("stats.txt", "w");
        fprintf(f_stats, "%d %d %d %d %fs", l, a, b, N, end - start);
        fclose(f_stats);
    }

    free(particles);

    MPI_Finalize();
    return 0;
}
