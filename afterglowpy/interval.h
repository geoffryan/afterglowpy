#ifndef AFTERGLOWPY_INTERVAL
#define AFTERGLOWPY_INTERVAL

/*
 * Here are 3 helper structs for performing global adaptive integration.
 * Each Mesh struct is a priority queue implemented with as a simple
 * binary heap. Elements may be added, and the element with the worst error can
 * be retrieved.
 *
 * Each mesh is a container of Intervals. An Interval is a primitive struct
 * contains only its left and right bounds, the value of the integral and its
 * error over that interval.  The Interval3 and Interval5 are identical to
 * Interval but contain 3 (or 5) extra slots to contain function memorized
 * function values.
 *
 * The Mesh3 and Mesh5 are priority queues of Interval3s and Interval5s
 * respectively.
 */

struct Interval
{
    double a;
    double b;
    double I;
    double err;
};
typedef struct Interval Interval;

struct Mesh
{
    size_t totalSize;
    size_t N;
    struct Interval *heap;
};
typedef struct Mesh Mesh;

struct Interval3
{
    double a;
    double b;
    double I;
    double err;
    double fa;
    double fb;
    double fm;
};
typedef struct Interval3 Interval3;

struct Mesh3
{
    size_t totalSize;
    size_t N;
    struct Interval3 *heap;
};
typedef struct Mesh3 Mesh3;

struct Interval5
{
    double a;
    double b;
    double I;
    double err;
    double fa;
    double fb;
    double fl;
    double fm;
    double fr;
};
typedef struct Interval5 Interval5;

struct Mesh5
{
    size_t totalSize;
    size_t N;
    struct Interval5 *heap;
};
typedef struct Mesh5 Mesh5;

void meshInit(struct Mesh *m);
void meshFree(struct Mesh *m);
void meshInsert(struct Mesh *m, struct Interval i);
struct Interval meshExtract(struct Mesh *m);
double meshTotalIntegral(struct Mesh *m);
double meshTotalError(struct Mesh *m);
void meshHeapifyUp(struct Mesh *m);
void meshHeapifyDown(struct Mesh *m);
int meshCheck(struct Mesh *m);

void mesh3Init(struct Mesh3 *m);
void mesh3Free(struct Mesh3 *m);
void mesh3Insert(struct Mesh3 *m, struct Interval3 i);
struct Interval3 mesh3Extract(struct Mesh3 *m);
double mesh3TotalIntegral(struct Mesh3 *m);
double mesh3TotalError(struct Mesh3 *m);
void mesh3HeapifyUp(struct Mesh3 *m);
void mesh3HeapifyDown(struct Mesh3 *m);
int mesh3Check(struct Mesh3 *m);

void mesh5Init(struct Mesh5 *m);
void mesh5Free(struct Mesh5 *m);
void mesh5Insert(struct Mesh5 *m, struct Interval5 i);
struct Interval5 mesh5Extract(struct Mesh5 *m);
double mesh5TotalIntegral(struct Mesh5 *m);
double mesh5TotalError(struct Mesh5 *m);
void mesh5HeapifyUp(struct Mesh5 *m);
void mesh5HeapifyDown(struct Mesh5 *m);
int mesh5Check(struct Mesh5 *m);

#endif
