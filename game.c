
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

#define DEAD 0
#define ALIVE 1

#define generations 14

#define matrix_size 20

/* do kazdego elementu macierzy stosujemy zasady gry i przechowujemy je w temp */
/* czyli tworzymy nowa generacje*/
void new_generation(int **matrix, int **temp,
                    int la, int lb)
{
    int i, j, sum;

    /* do kazdego elementu macierzy stosujemy zasady gry i przechowujemy je w temp */
    for (i = la; i <= lb; i++)
    {
        for (j = 1; j <= matrix_size; j++)
        {

            /* szukamy sasiadow danej komorki macierzy */
            sum = matrix[i - 1][j - 1] + matrix[i - 1][j] + matrix[i - 1][j + 1] + matrix[i][j - 1] + matrix[i][j + 1] + matrix[i + 1][j - 1] + matrix[i + 1][j] + matrix[i + 1][j + 1];

            /* sprawdzamy czy tworzy sie nowa lub jakas umiera */

            /* umiera z izolacji lub bo ma za duzo sasiadow */
            if (sum < 2 || sum > 3)
                temp[i][j] = DEAD;

            /* tworzy sie nowa */
            else if (sum == 3)
                temp[i][j] = ALIVE;

            /* w pozostalych wypadkach zpstawiamy */
            else
                temp[i][j] = matrix[i][j];
        }
    }
}

/* podmieniamy maciez temp na matrix*/
void swap_matrices(int **matrix, int **temp,
                   int *changes, int la, int lb)
{
    int i, j, aux;

    /* matrix = temp */
    *changes = 0;
    for (i = la; i <= lb; i++)
    {
        for (j = 1; j <= matrix_size; j++)
        {
            aux = matrix[i][j];
            matrix[i][j] = temp[i][j];
            temp[i][j] = aux;
            if (matrix[i][j] != temp[i][j])
                (*changes)++;
        }
    }
}

/* tutaj wypisujemy nasza maciez */
void output_life_matrix(int **matrix,
                        int iter, int changes)
{
    int i, j;

    /* klasyczne wypisywanie 2d macierzy */
    for (i = 0; i < matrix_size; i++)
    {
        for (j = 0; j < matrix_size; j++)
        {
            printf("%d ", matrix[i][j]);
        }

        printf("\n");
    }
    fprintf(stderr, "Liczba iteracji: %d  ---- Liczba zmian: %d\n",
            iter, changes);
}

/* rownomiernie rozprzestrzeniamy random */
double rand01()
{
    return (double)rand() / (double)RAND_MAX;
}

/* inicjujemy macierz o wyznaczonej wielkosci*/
void initial(int **matrix, int **temp, int myrank)
{
    int i, j;

    /* zapelniamy 0 */
    for (j = 0; j < matrix_size + 2; j++)
    {
        matrix[0][j] = DEAD;
        matrix[matrix_size + 1][j] = DEAD;
        temp[0][j] = DEAD;
        temp[matrix_size + 1][j] = DEAD;
    }
    for (i = 0; i < matrix_size + 2; i++)
    {
        matrix[i][0] = DEAD;
        matrix[i][matrix_size + 1] = DEAD;
        temp[i][0] = DEAD;
        temp[i][matrix_size + 1] = DEAD;
    }
    if (myrank == 0)
        fprintf(stderr, "wyznaczona\n");

    /* inicjalizujemy */
    for (i = 1; i <= matrix_size; i++)
    {
        for (j = 1; j <= matrix_size; j++)
            if (rand01() > 0.5)
                matrix[i][j] = ALIVE;
            else
                matrix[i][j] = DEAD;
    }
    if (myrank == 0)
        fprintf(stderr, "Population initialized\n");
}

/* dzielimy na strips, ktore bazuja na wierszach */
void decompose_domain(int *start_strip, int *end_strip,
                      int MPIsize, int myrank)
{
    int ipointer, pe, strip_size;

    if (myrank == 0)
        fprintf(stderr, "Strips: ");
    start_strip[1] = 1;
    strip_size = (matrix_size) / (MPIsize - 1) + 1;
    for (pe = 1; pe < MPIsize; pe++)
    {
        end_strip[pe] = start_strip[pe] + strip_size - 1;
        if (pe == MPIsize - 1)
            end_strip[pe] = matrix_size;
        else
            start_strip[pe + 1] = end_strip[pe] + 1;
        if (myrank == 0)
        {
            fprintf(stderr, " [ %d , %d ]",
                    start_strip[pe], end_strip[pe]);
            if (pe == MPIsize - 1)
                fprintf(stderr, "\n");
        }
    }
}

void strips_boundary_exchange(int **matrix, int *start_strip,
                              int *end_strip, int myrank,
                              int MPIsize)
{
    MPI_Status recv_status;
    int top_process, bot_process;

    top_process = myrank + 1;
    bot_process = myrank - 1;

    /* pracujemy na nieparzystych */
    if (myrank % 2 != 0)
    {
        if (myrank < MPIsize - 1) /* wysylamy do prawo- gora */
            MPI_Send(&matrix[end_strip[myrank]][0], matrix_size + 2,
                     MPI_INT, top_process, 102, MPI_COMM_WORLD);
        if (myrank > 1) /* otrzymujemy z lewo dol */
            MPI_Recv(&matrix[start_strip[myrank] - 1][0], matrix_size + 2,
                     MPI_INT, bot_process, 101, MPI_COMM_WORLD,
                     &recv_status);
        if (myrank > 1) /* wysylamy do lewo dol to */
            MPI_Send(&matrix[start_strip[myrank]][0], matrix_size + 2,
                     MPI_INT, bot_process, 104, MPI_COMM_WORLD);
        if (myrank < MPIsize - 1) /* otrzymujemy z prawo- gora */
            MPI_Recv(&matrix[end_strip[myrank] + 1][0], matrix_size + 2,
                     MPI_INT, top_process, 103, MPI_COMM_WORLD,
                     &recv_status);
    }
    /* parzyste */
    else
    {
        if (myrank > 1) /* otrzymujemy z lewo dol */
            MPI_Recv(&matrix[start_strip[myrank] - 1][0], matrix_size + 2,
                     MPI_INT, bot_process, 102, MPI_COMM_WORLD,
                     &recv_status);
        if (myrank < MPIsize - 1) /* wysylamy do prawo- gora */
            MPI_Send(&matrix[end_strip[myrank]][0], matrix_size + 2,
                     MPI_INT, top_process, 101, MPI_COMM_WORLD);
        if (myrank < MPIsize - 1) /* otrzymujemy z prawo- gora */
            MPI_Recv(&matrix[end_strip[myrank] + 1][0], matrix_size + 2,
                     MPI_INT, top_process, 104, MPI_COMM_WORLD,
                     &recv_status);
        if (myrank > 1) /* wysylamy do lewo dol to */
            MPI_Send(&matrix[start_strip[myrank]][0], matrix_size + 2,
                     MPI_INT, bot_process, 103, MPI_COMM_WORLD);
    }
}

int main(int argc, char *argv[])
{
    int i, j, iter, changes, total_changes;
    int **matrix, **temp;

    int *transfer;
    int *start_strip, *end_strip;
    int myrank, MPIsize, pe;
    int ipointer, strip_size, mesg_size, from_strip;

    MPI_Status recv_status;

    fprintf(stderr, "Game starts\n");

    /* init MPI */
    MPI_Init(&argc, &argv);
    /* mmyrank -  id */
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    /* ile mozemy miec procesow*/
    MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);

    fprintf(stderr, "MPI communicator has been attached\n");

    /* alokujemy pamiec dla macierzy */
    matrix = (int **)malloc(sizeof(int *) * (matrix_size + 2));
    temp = (int **)malloc(sizeof(int *) * (matrix_size + 2));
    for (i = 0; i < matrix_size + 2; i++)
    {
        matrix[i] = (int *)malloc(sizeof(int) * (matrix_size + 2));
        temp[i] = (int *)malloc(sizeof(int) * (matrix_size + 2));
    }
    transfer = (int *)malloc(sizeof(int) *
                             (matrix_size + 2) * (matrix_size + 2));
    if (myrank == 0)
        fprintf(stderr, "Memory for matrix set aside \n");

    /* ustalamy nasze strips (wielkosc)*/
    start_strip = (int *)malloc(sizeof(int) * MPIsize);
    end_strip = (int *)malloc(sizeof(int) * MPIsize);
    decompose_domain(start_strip, end_strip, MPIsize, myrank);

    /* inicjalizujemy granice i macierze         */
    /* node 0 jest tylko jako output*/
    initial(matrix, temp, myrank);

    /* pierwsza macierz */
    if (myrank == 0)
    {
        printf("\n");
        output_life_matrix(matrix, 0, 0);
    }

    /* iterujemy po generacjach */
    for (iter = 1; iter <= generations; iter++)
    {
        /* nie moga byc 0 */
        if (myrank != 0)
        {
            /* generujemy nowa generavje w strip */
            new_generation(matrix, temp,
                           start_strip[myrank],
                           end_strip[myrank]);

            /* podmieniamy stara macierz na nowa*/
            swap_matrices(matrix, temp, &changes,
                          start_strip[myrank],
                          end_strip[myrank]);

            /* management granic miedzy stripami */
            /* prawdziwe nie sa tu zmieniane */
            strips_boundary_exchange(matrix, start_strip, end_strip,
                                     myrank, MPIsize);

            /* synchronizacja procesow */
            MPI_Barrier(MPI_COMM_WORLD);

            /* liczymy liczbe zmian */
            MPI_Reduce(&changes, &total_changes, 1, MPI_INT,
                       MPI_SUM, 0, MPI_COMM_WORLD);

            /* wysylamy strip do node 0 */
            ipointer = 0;
            for (i = start_strip[myrank]; i <= end_strip[myrank]; i++)
                for (j = 0; j <= matrix_size + 1; j++)
                {
                    transfer[ipointer] = matrix[i][j];
                    ipointer++;
                }
            mesg_size = ipointer++;
            MPI_Send(&transfer[0], mesg_size,
                     MPI_INT, 0, 121, MPI_COMM_WORLD);
        }

        /* node 0 */
        else
        {
            /* synchronizacja wszystkich nodow */
            MPI_Barrier(MPI_COMM_WORLD);

            /* liczymy zmiany */
            MPI_Reduce(&changes, &total_changes, 1, MPI_INT,
                       MPI_SUM, 0, MPI_COMM_WORLD);

            /* otrzymujemy stripy */
            strip_size = (matrix_size) / (MPIsize - 1) + 1;
            mesg_size = (matrix_size + 2) * strip_size;
            for (pe = 1; pe <= MPIsize - 1; pe++)
            {
                MPI_Recv(&transfer[0], mesg_size,
                         MPI_INT, MPI_ANY_SOURCE, 121, MPI_COMM_WORLD,
                         &recv_status);
                from_strip = recv_status.MPI_SOURCE;
                ipointer = 0;
                for (i = start_strip[from_strip];
                     i <= end_strip[from_strip]; i++)
                {
                    for (j = 0; j <= matrix_size + 1; j++)
                    {
                        matrix[i][j] = transfer[ipointer];
                        ipointer++;
                    }
                }
            }

            /* wyswietlamy nowa generacje*/
            printf("\n");
            output_life_matrix(matrix, iter, total_changes);
        }
    }

    MPI_Finalize();
    exit(0);
}