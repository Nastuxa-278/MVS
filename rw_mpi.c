#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <pthread.h>
#include <assert.h>

#define L 0
#define R 1
#define U 2
#define D 3
#define OK 4

typedef unsigned int uint;

typedef struct particle {

    int x;
    int y;
    int rest_count;
    uint seed;
} particle;

typedef struct arguments{
  int rang;
  int count_proc;
  int l;
  int a;
  int b;
  int n;
  int N;
  double pl;
  double pr;
  double pu;
  double pd;
}arguments;

void walk(int rang, int count_proc, int l, int a, int b, int n, int N, double pl, double pr, double pu, double pd);

int choose_direction(double pl, double pr, double pu, double pd, uint* seed){
    double number = ((double) rand_r(seed));


    //divide length into 4 parts and choose destination
    if (number <= RAND_MAX*pl){
        return L;
    }

    else if (number <= RAND_MAX*(pl+pr)){
        return R;
    }

    else if(number <= RAND_MAX*(pl+pr+pu)){
        return U;
    }

    else{
        return D;
    }
}

void step(double pl, double pr, double pu, double pd, particle* part){
    int direction = choose_direction(pl, pr, pu, pd, &(part->seed));

    if(direction == L){
        part->x--;
    }
    else if (direction == R) {
        part->x++;
    }
    else if (direction == U) {
        part->y--;
    }
    else {
        part->y++;
    }
    part->rest_count--;
}


int check_bounds(particle* part, int l) {
    if(part->x < 0) {
        part->x = l - 1;
        return L;
    }
    if(part->x >= l) {
        part->x = 0;
        return R;
    }
    if (part->y < 0) {
        part->y = l - 1;
        return U;
    }
    if (part->y >= l) {
        part->y = 0;
        return D;
    }
    return OK;
}


void insert(particle* particles, int* particles_count, int* particles_real_len, particle* from, int* from_len){
    for(int i = 0; i < (*from_len); ++i) {
        if (*particles_count >= *particles_real_len) {
            (*particles_real_len) = 2*(*particles_real_len);
            particles = (particle*) realloc(particles, (*particles_real_len) * sizeof(particle));
        }
        particles[*particles_count] = from[i];
        (*particles_count) = (*particles_count) + 1;
    }
}


void neighbour_rangs(int rang, int count_proc, int l, int a, int b, int x, int y, int* nb_rangs){
    int rang_l, rang_r, rang_u, rang_d;

    if (x == 0) {
        rang_l = rang + a - 1;
    }
    else {
        rang_l = rang - 1;
    }
    if (x == a - 1) {
        rang_r = rang - a + 1;
    }
    else{
        rang_r = rang + 1;
    }
    if (y == 0) {
        rang_u = count_proc - a + rang;
    }
    else {
        rang_u = rang - a;
    }
    if(y == b - 1){
        rang_d = rang + a - count_proc;
    }
    else{
        rang_d = rang + a;
    }
    nb_rangs[0] = rang_l;
    nb_rangs[1] = rang_r;
    nb_rangs[2] = rang_u;
    nb_rangs[3] = rang_d;
}

void* walk_void_args(void* args_){
  arguments* args = (arguments*) args_;
  walk(args->rang, args->count_proc, args->l, args->a, args->b, args->n, args->N, args->pl, args->pr, args->pu, args->pd);

  return 0;

}


void walk(int rang, int count_proc, int l, int a, int b, int n, int N, double pl, double pr, double pu, double pd) {
    double start = MPI_Wtime();


    int rang_x = rang % a;
    int rang_y = rang / a;

    int nb_rangs[4];
    neighbour_rangs(rang, count_proc, l, a, b, rang_x, rang_y, nb_rangs);

    int particles_real_len = 2*N;
    int particles_count = N;
    particle* particles = (particle*) malloc(sizeof(particle) * particles_real_len);


    uint seed;
    uint* seeds = (uint*)malloc(sizeof(uint)* count_proc);
    if (rang == 0) {
      srand(time(NULL));
      for (int i = 0; i < count_proc; ++i) {
          seeds[i] = rand();
      }
    }
    MPI_Scatter(seeds, 1, MPI_UNSIGNED, &seed, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    free(seeds);


    for(int i = 0; i < N; ++i) {
        particles[i].x = rand_r(&seed) % l;
        particles[i].y = rand_r(&seed) % l;
        particles[i].rest_count = n;
        particles[i].seed = rand_r(&seed);
    }


    particle* to_left = (particle*) malloc(sizeof(particle) * N);
    particle* to_right = (particle*) malloc(sizeof(particle) * N);
    particle* to_up = (particle*) malloc(sizeof(particle) * N);
    particle* to_down = (particle*) malloc(sizeof(particle) * N);
    int to_left_len = 0;
    int to_right_len = 0;
    int to_up_len = 0;
    int to_down_len = 0;

    int to_left_real_len = N;
    int to_right_real_len = N;
    int to_up_real_len = N;
    int to_down_real_len = N;


    int count_finished = 0;


    int count_iteration = 20;
    int is_end = 0;

    while (!is_end){
        to_left_len = 0;
        to_right_len = 0;
        to_up_len = 0;
        to_down_len = 0;
        for(int k = 0; k < count_iteration; ++k){


            int i = 0;
            while (i < particles_count){


                if (particles[i].rest_count == 0) {
                    count_finished++;
                    particles[i] = particles[particles_count - 1];
                    particles_count--;
                    i--;
                }

                else {

                    step(pl, pr, pu, pd, particles + i);


                    int res = check_bounds(particles + i, l);

                    if (res==L){
                        if(to_left_len >= to_left_real_len) {
                            to_left_real_len = 2* to_left_real_len;
                            to_left = (particle*)realloc(to_left, to_left_real_len*sizeof(particle));
                        }
                        to_left[to_left_len] = particles[i];
                        particles[i] = particles[particles_count - 1];
                        particles_count--;
                        to_left_len++;
                        i--;
                    }
                    else if (res==R){
                        if(to_right_len >= to_right_real_len) {
                            to_right_real_len = 2* to_right_real_len;
                            to_right = (particle*)realloc(to_right, to_right_real_len*sizeof(particle));
                        }
                        to_right[to_right_len] = particles[i];
                        particles[i] = particles[particles_count - 1];
                        particles_count--;
                        to_right_len++;
                        i--;
                    }
                    else if (res==U){
                        if(to_up_len >= to_up_real_len) {
                            to_up_real_len = 2* to_up_real_len;
                            to_up = (particle*)realloc(to_up, to_up_real_len*sizeof(particle));
                        }
                        to_up[to_up_len] = particles[i];
                        particles[i] = particles[particles_count - 1];
                        particles_count--;
                        to_up_len++;
                        i--;
                    }
                    else if (res==D){
                        if(to_down_len >= to_down_real_len) {
                            to_down_real_len = 2* to_down_real_len;
                            to_down= (particle*)realloc(to_down, to_down_real_len*sizeof(particle));
                        }
                        to_down[to_down_len] = particles[i];
                        particles[i] = particles[particles_count - 1];
                        particles_count--;
                        to_down_len++;
                        i--;
                    }
                }
                i++;
            }
        }



        MPI_Request requests_to[4];
        MPI_Request requests_from[4];


        int from_left_len = 0;
        int from_right_len = 0;
        int from_up_len = 0;
        int from_down_len = 0;


        //use both MPI_Isend and MPI_Irecv
        MPI_Issend(&to_left_len, 1, MPI_INT, nb_rangs[0], L, MPI_COMM_WORLD, requests_from);
        MPI_Irecv(&from_left_len, 1, MPI_INT, nb_rangs[0], R, MPI_COMM_WORLD, requests_to);

        MPI_Issend(&to_right_len, 1, MPI_INT, nb_rangs[1], R, MPI_COMM_WORLD, requests_from + 1);
        MPI_Irecv(&from_right_len, 1, MPI_INT, nb_rangs[1], L, MPI_COMM_WORLD, requests_to+1);

        MPI_Issend(&to_up_len, 1, MPI_INT, nb_rangs[2], U, MPI_COMM_WORLD, requests_from + 2);
        MPI_Irecv(&from_up_len, 1, MPI_INT, nb_rangs[2], D, MPI_COMM_WORLD, requests_to+2);

        MPI_Issend(&to_down_len, 1, MPI_INT, nb_rangs[3], D, MPI_COMM_WORLD, requests_from + 3);
        MPI_Irecv(&from_down_len, 1, MPI_INT, nb_rangs[3], U, MPI_COMM_WORLD, requests_to+3);


        MPI_Status status_to[4];
        MPI_Status status_from[4];


        for(int j = 0; j <4; ++j) {
            MPI_Wait(requests_from + j, status_to + j);
            MPI_Wait(requests_to + j, status_from + j);
        }


        particle* from_left = (particle*) malloc(sizeof(particle) * from_left_len);
        particle* from_right = (particle*) malloc(sizeof(particle) * from_right_len);
        particle* from_up = (particle*) malloc(sizeof(particle) * from_up_len);
        particle* from_down = (particle*) malloc(sizeof(particle) * from_down_len);

        MPI_Request reqfrom[4];
        MPI_Request reqto[4];
        MPI_Status statto[4];
        MPI_Status statfrom[4];



        MPI_Isend(to_left, to_left_len*sizeof(particle),MPI_BYTE, nb_rangs[0], L, MPI_COMM_WORLD, reqfrom);
        MPI_Irecv(from_left, from_left_len*sizeof(particle), MPI_BYTE, nb_rangs[0], R, MPI_COMM_WORLD, reqto);

        MPI_Isend(to_right, to_right_len*sizeof(particle), MPI_BYTE, nb_rangs[1], R, MPI_COMM_WORLD, reqfrom+1);
        MPI_Irecv(from_right, from_right_len*sizeof(particle), MPI_BYTE, nb_rangs[1], L, MPI_COMM_WORLD, reqto+1);

        MPI_Isend(to_up, to_up_len*sizeof(particle), MPI_BYTE, nb_rangs[2], U, MPI_COMM_WORLD, reqfrom+2);
        MPI_Irecv(from_up, from_up_len*sizeof(particle), MPI_BYTE, nb_rangs[2], D, MPI_COMM_WORLD, reqto+2);

        MPI_Isend(to_down, to_down_len*sizeof(particle),MPI_BYTE, nb_rangs[3], D, MPI_COMM_WORLD, reqfrom+3);
        MPI_Irecv(from_down, from_down_len*sizeof(particle), MPI_BYTE, nb_rangs[3], U, MPI_COMM_WORLD, reqto+3);

        for(int j = 0; j <4; ++j) {
            MPI_Wait(reqfrom + j, statto + j);
            MPI_Wait(reqto + j,  statfrom + j);
        }


        insert(particles, &particles_count, &particles_real_len, from_left, &from_left_len);
        insert(particles, &particles_count, &particles_real_len, from_right, &from_right_len);
        insert(particles, &particles_count, &particles_real_len, from_up, &from_up_len);
        insert(particles, &particles_count, &particles_real_len, from_down, &from_down_len);



        //use MPI_Reduce to save in all_count_finished whole sum of finished particles
        int all_count_finished = 0;
        MPI_Reduce(&count_finished, &all_count_finished, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);


        MPI_Request* requests_send = (MPI_Request*)malloc(sizeof(MPI_Request) * count_proc);
        MPI_Request* requests_recv = (MPI_Request*)malloc(sizeof(MPI_Request) * count_proc);


        if (rang == 0) {
          if(all_count_finished == count_proc * N){
            is_end = 1;
          }
          else{
            is_end = 0;
          }

          for(int j = 0; j < count_proc; ++j) {
              MPI_Isend(&is_end, 1 ,MPI_INT, j, 0, MPI_COMM_WORLD, requests_send + j);
          }
        }

        MPI_Irecv(&is_end, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, requests_recv + rang);


        MPI_Barrier(MPI_COMM_WORLD);


        free(requests_send);
        free(requests_recv);

        free(from_left);
        free(from_right);
        free(from_up);
        free(from_down);

    }

    double end = MPI_Wtime();

    int* result = (int*)malloc(count_proc * 2 * sizeof(int));


    int rang_and_count[2];
    rang_and_count[0] = rang;
    rang_and_count[1] = count_finished;


    MPI_Gather(rang_and_count, 2, MPI_INT, result, 2, MPI_INT, 0, MPI_COMM_WORLD);


    if(rang==0) {
      double program_time = end - start;
      FILE *f_stats;
      f_stats = fopen("stats.txt", "w");
      fprintf(f_stats, "%d %d %d %d %d %f %f %f %f %f\n", l, a, b, n, N, pl, pr, pu, pd, program_time);

      for (int i = 0; i < 2*count_proc; i = i + 2){
        fprintf(f_stats, "%d : %d\n", result[i], result[i+1]);
      }
      fclose(f_stats);
    }

    free(result);
    free(particles);
    free(to_left);
    free(to_right);
    free(to_up);
    free(to_down);
}


int main(int argc, char * argv[]) {
    MPI_Init(&argc, &argv);

    int rang, count_proc;

    MPI_Comm_rank(MPI_COMM_WORLD, &rang);
    MPI_Comm_size(MPI_COMM_WORLD, &count_proc);

    int l, a, b, n, N;
    double pl, pr, pu, pd;

    l = atoi(argv[1]);
    a = atoi(argv[2]);
    b = atoi(argv[3]);
    n = atoi(argv[4]);
    N = atoi(argv[5]);
    pl = atof(argv[6]);
    pr = atof(argv[7]);
    pu = atof(argv[8]);
    pd = atof(argv[9]);


    assert (count_proc == a * b);

    arguments args;
    args.rang = rang;
    args.count_proc = count_proc;
    args.l = l;
    args.a = a;
    args.b = b;
    args.n = n;
    args.N = N;
    args.pl = pl;
    args.pr = pr;
    args.pu = pu;
    args.pd = pd;

    walk_void_args(&args);

    MPI_Finalize();
    return 0;
}
