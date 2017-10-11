#include "particles.h"
#include "myrand.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define PRINT_OPAQUE_STRUCT(p)  print_mem((p), sizeof(*(p)))

void print_mem(void const *vp, size_t n)
{
    unsigned char const *p = vp;
    for (size_t i=0; i<n; i++)
        printf("%02x%s", p[i], i==0 ? "" : i % 16 == 15 ? "\n" : i % 4 == 3 ? " " : "");
    putchar('\n');
};

//ldoc
/**
 * # Particle data stucture implementation
 *
 * The particle data structure is one of the places where you may
 * want to improve the implementation.  There are several aspects
 * that are not entirely satisfactory; are we better off with the
 * parallel array layout we have now (which is good for vectorization,
 * but perhaps not locality), or with an array-of-structs layout?
 * Should the particle type (which is really only a single bit of
 * information) be stored as a full integer, or as some smaller data
 * type?  Is there a better way of handling the bin data structures
 * (noting that this is not the bin data structure from the original
 * code)?  Play with things and find out!
 * 
 * ## Allocate and free particle data
 */
/*
particles_t* alloc_particles_t(int nbinx, int N, float L)
{
    int size_total = nbinx*nbinx;
    particles_t* p = malloc(sizeof(particles_t));
    p->nbinx = nbinx;
    p->N = N;
    p->L = L;
    p->type =  (int*)   calloc(  N, sizeof(int));
    p->x =     (float*) calloc(2*N, sizeof(float));
    p->v =     (float*) calloc(2*N, sizeof(float));
    p->f =     (float*) calloc(2*N, sizeof(float));
    p->next =  (int*)   calloc(  N, sizeof(int));
    p->cells = (int*)   calloc(size_total, sizeof(int));
    return p;
}
*/

/*
void free_particles_t(particles_t* p)
{
    free(p->cells);
    free(p->next);
    free(p->f);
    free(p->v);
    free(p->type);
    free(p);
}
*/

/**
 * ## Neighbor list computations
 *
 * Right now, computing a neighbor list means computing the `cells`
 * and `next` arrays; we don't touch the position, velocity, or force
 * arrays.  But perhaps we should!  As currently written, traversing
 * a neighbor list involves jumping all over memory to grab the information
 * for different particles; you may want to consider something that makes
 * better use of cache locality.
 */
/*
void compute_nbr_lists(particles_t* particles)
{
    int* restrict cells = particles->cells;
    int* restrict next  = particles->next;
    float* restrict x = particles->x;
    float L = particles->L;
    int N = particles->N;
    int nbinx = particles->nbinx;

    // Recompute neighbor list
    const int size_total = nbinx*nbinx;
    memset(cells, -1, size_total * sizeof(int));
    for (int i=0; i<N; i++) {
        int binx = coord_to_index(x[2*i+0], nbinx, L);
        int biny = coord_to_index(x[2*i+1], nbinx, L);
        int t = binx + biny*nbinx;
        next[i] = cells[t];
        cells[t] = i;
    }
}
*/

/**
 * ## Initialization
 *
 * The `initial_circle` function puts all the active moshers in a big
 * circle in the middle of the domain at the start; the `init_ric`
 * function distributes the active moshers uniformly through the
 * domain and assigns them a random velocity as well as position.
 * It might be interesting to see how the behavior depends on this
 * initial distribution.
 */

/*
void init_ric(particles_t* particles, float speed)
{
    float* restrict x = particles->x;
    float* restrict v = particles->v;
    int* restrict type = particles->type;
    int N = particles->N;
    float L = particles->L;

    for (int i=0; i<N; i++) {
        float t = 2*M_PI*ran_ran2();

        x[2*i+0] = L*ran_ran2();
        x[2*i+1] = L*ran_ran2();

        if (ran_ran2() > 0.16){
            v[2*i+0] = 0.0;
            v[2*i+1] = 0.0;
            type[i] = BLACK;
        } else {
            v[2*i+0] = speed * sin(t);
            v[2*i+1] = speed * cos(t);
            type[i] = RED;
        }
    }
}


void init_circle(particles_t* particles, float speed)
{
    float* restrict x = particles->x;
    float* restrict v = particles->v;
    int* restrict type = particles->type;
    int N = particles->N;
    float L = particles->L;
    
    for (int i=0; i<N; i++){
        float tx = L*ran_ran2();
        float ty = L*ran_ran2();
        float tt = 2*M_PI*ran_ran2();

        x[2*i+0] = tx;
        x[2*i+1] = ty;

        // the radius for which 30% of the particles are red on avg
        float dd2 = (tx-L/2)*(tx-L/2) + (ty-L/2)*(ty-L/2);
        float rad2 = 0.16*L*L / M_PI;

        if (dd2 < rad2)
            type[i] = RED;
        else
            type[i] = BLACK;

        if (type[i] == RED) {
            v[2*i+0] = speed*cos(tt);
            v[2*i+1] = speed*sin(tt);
        } else {
            v[2*i+0] = 0.0;
            v[2*i+1] = 0.0;
        }
    }
}
*/


/* CELLS */

particles_t* alloc_particles_t(int nbinx, int N, float L)
{
    int size_total = nbinx*nbinx;
    particles_t* p = malloc(sizeof(particles_t));
    p->nbinx = nbinx;
    p->N = N;
    p->L = L;
    p->particles =  (particle_t*)   calloc(  N, sizeof(particle_t));

    particle_t** cells = calloc(size_total, sizeof(particle_t));
    //p->cells = (particle_t*)   calloc(size_total, sizeof(particle_t));
    for(int i=0;i<size_total;i++)
    {
        cells[i] = 0;
    }
    p->cells=cells;
    return p;
}

void free_particles_t(particles_t* p)
{
    free(p->cells);
    free(p->particles);
    free(p);
}


void compute_nbr_lists(particles_t* particles)
{
    //printf("---Compute Neighbor Lists\n");
    //fflush(stdout);

    int N = particles->N;
    int nbinx = particles->nbinx;
    float L = particles->L;
    particle_t *p =  malloc(sizeof(particle_t));

    // Recompute neighbor list
    const int size_total = nbinx*nbinx;
    memset(particles->cells, -1, size_total * sizeof(particle_t));
    for(int i=0;i<N;i++) 
    {        
        p = &particles->particles[i];

        int binx = coord_to_index(p->x, nbinx, L);
        int biny = coord_to_index(p->y, nbinx, L);
        int t = binx + biny*nbinx;
        
        //printf("Particle [%i] (%f,%f) Bin [%i]\n",i,p->x,p->y,t);
        //fflush(stdout);

        p->next = particles->cells[t];
        particles->cells[t] = p;        
        
        //printf("Assigned %p to %p\n",p->next,particles->cells[t]);
        //fflush(stdout);
    }
}

/*
void compute_nbr_lists(particles_t* particles)
{
    int* restrict cells = particles->cells;
    int* restrict next  = particles->next;
    float* restrict x = particles->x;
    float L = particles->L;
    int N = particles->N;
    int nbinx = particles->nbinx;

    // Recompute neighbor list
    const int size_total = nbinx*nbinx;
    memset(cells, -1, size_total * sizeof(int));
    for (int i=0; i<N; i++) {
        int binx = coord_to_index(x[2*i+0], nbinx, L);
        int biny = coord_to_index(x[2*i+1], nbinx, L);
        int t = binx + biny*nbinx;
        next[i] = cells[t];
        cells[t] = i;
    }
}
*/

void init_ric(particles_t* particles, float speed)
{
    printf("Init ric\n");
    fflush(stdout);
    
    float x;
    float y;
    float vx;
    float vy;
    int type;

    int binx;
    int biny;
    int t;

    int N = particles->N;
    int nbinx = particles->nbinx;
    float L = particles->L;


    for (int i=0; i<N; i++) {
        float t = 2*M_PI*ran_ran2();

        x = L*ran_ran2();
        y = L*ran_ran2();

        if (ran_ran2() > 0.16){
            vx = 0.0;
            vy = 0.0;
            type = BLACK;
        } else {
            vx = speed * sin(t);
            vy = speed * cos(t);
            type = RED;
        }

        particle_t *p =  malloc(sizeof(particle_t));
        p->x = x;
        p->y = y;
        p->vx = vx;
        p->vy = vy;
        p->fy = 0.0;
        p->fy = 0.0;
        p->type=type;
        p->next = NULL;

        particles->particles[i] = *p;


        /* Needed?
        binx = coord_to_index(x, nbinx, L);
        biny = coord_to_index(y, nbinx, L);
        t = binx + biny*nbinx;
        p->next = p->cells[(int)t];
        p->cells[(int)t] = p;
        */
        //printf("X=%f Y=%f at (%i,%i) in bin %i \n",p->x,p->y,binx,biny,(int)t);
        //fflush(stdout);
    }
}

void init_circle(particles_t* particles, float speed)
{
    int N = particles->N;
    float L = particles->L;

    float x;
    float y;
    float vx;
    float vy;
    int type;

    int binx;
    int biny;
    int t;
    
    for (int i=0; i<N; i++){
        float x = L*ran_ran2();
        float y = L*ran_ran2();
        float tt = 2*M_PI*ran_ran2();


        // the radius for which 30% of the particles are red on avg
        float dd2 = (x-L/2)*(x-L/2) + (y-L/2)*(y-L/2);
        float rad2 = 0.16*L*L / M_PI;

        if (dd2 < rad2)
            type = RED;
        else
            type = BLACK;

        if (type == RED) {
            vx = speed*cos(tt);
            vy = speed*sin(tt);
        } else {
            vx = 0.0;
            vy = 0.0;
        }

        particle_t *p = (particle_t*)  calloc(1, sizeof(particle_t));
        p->x = x;
        p->y = y;
        p->vx = vx;
        p->vy = vy;
        p->fy = 0.0;
        p->fy = 0.0;
        p->type=type;
        p->next = NULL;

        particles->particles[i] = *p;

        /* Needed?
        binx = coord_to_index(x, nbinx, L);
        biny = coord_to_index(y, nbinx, L);
        t = binx + biny*nbinx;
        p->next = p->cells[(int)t];
        p->cells[(int)t] = p;
        */
        //printf("X=%f Y=%f at (%i,%i) in bin %i \n",p->x,p->y,binx,biny,(int)t);
        //fflush(stdout);
    }
}


FILE* start_frames(const char* fname)
{
    FILE* fp = fopen(fname, "w");
    if (fp == NULL) {
        fprintf(stderr, "Could not open %s for output\n", fname);
        exit(-1);
    }
    fprintf(fp, "PTag,PId,PLocX,PLocY,PDirX,PDirY\n");
    return fp;
}

void end_frames(FILE* fp)
{
    fclose(fp);
}

void write_frame(FILE* fp, particles_t* particles)
{
    int n = particles->N;
    float L = particles->L;
    for (int i = 0; i < n; ++i) {
        fprintf(fp, "%d,%d,%g,%g,%g,%g\n", particles->particles[i].type, i+1,
                (particles->particles[i].x)/L, (particles->particles[i].y)/L,
                (particles->particles[i].vx)/L, (particles->particles[i].vy)/L);
    }
}