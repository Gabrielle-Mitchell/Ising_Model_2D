#include <stdbool.h>
typedef struct {int x,y;} lat_type;
void initialize(int size, int lat[size+1][size+1]);
void output(int size, int lat[size+1][size+1]);
void choose_random_pos_lat(lat_type *);
int energy_pos(lat_type,int size, int lat[size+1][size+1]);
bool test_flip(lat_type pos, int *de,int size, int lat[size+1][size+1]);
void flip(lat_type pos, int size, int lat[size+1][size+1]);
void transient_results(int size, int lat[size+1][size+1]);
int total_magnetization(int size, int lat[size+1][size+1]);
int total energy(int size, int lat[size+1][size+1]);
int find_max(int nintervals, float C[]);