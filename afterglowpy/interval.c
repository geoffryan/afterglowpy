#include <stdio.h>
#include <stdlib.h>
#include "interval.h"

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

/***  Mesh  ***/

void meshInit(struct Mesh *m)
{
    size_t size = 4;
    m->totalSize = size;
    m->N = 0;
    m->heap = (struct Interval *)malloc(size * sizeof(struct Interval));
}

void meshFree(struct Mesh *m)
{
    m->totalSize = 0;
    m->N = 0;
    free(m->heap);
    m->heap = NULL;
}

void meshInsert(struct Mesh *m, struct Interval *i)
{
    // Resize if necessary
    while(m->N >= m->totalSize)
    {
        m->totalSize *= 2;
        m->heap = (struct Interval *)realloc(m->heap,
                                    m->totalSize * sizeof(struct Interval));
    }
    // Add interval to end of heap
    m->heap[m->N] = *i;
    (m->N)++;

    // Restore ordering
    meshHeapifyUp(m);
}

void meshExtract(struct Mesh *m, struct Interval *worst)
{
    *worst = m->heap[0];

    m->heap[0] = m->heap[m->N-1];
    (m->N)--;

    meshHeapifyDown(m);

}

double meshTotalIntegral(struct Mesh *m)
{
    double I = 0.0;

    size_t i = 0;
    for(i=0; i<m->N; i++)
        I += m->heap[i].I;
    
    return I;
}

double meshTotalError(struct Mesh *m)
{
    double err = 0.0;

    size_t i = 0;
    for(i=0; i<m->N; i++)
        err += m->heap[i].err;
    
    return err;
}

void meshHeapifyUp(struct Mesh *m)
{
    size_t c = m->N-1;
    size_t p = (c-1)/2;

    while(c != 0 && m->heap[p].err < m->heap[c].err)
    {
        struct Interval tempP = m->heap[p];
        m->heap[p] = m->heap[c];
        m->heap[c] = tempP;
        c = p;
        p = (c-1)/2;
    }
}

void meshHeapifyDown(struct Mesh *m)
{
    size_t p = 0;
    size_t c1 = 2*p + 1;
    size_t c2 = c1 + 1;


    while(c1 < m->N)
    {
        //Find the child with largest error
        size_t c = c1;
        double e = m->heap[c1].err;

        if(c2 < m->N && m->heap[c2].err > e)
        {
            c = c2;
            e = m->heap[c2].err;
        }

        // If the child is already in order then we're done.
        if(e <= m->heap[p].err)
            break;

        // Otherwise, swap with the child.
        struct Interval tempP = m->heap[p];
        m->heap[p] = m->heap[c];
        m->heap[c] = tempP;
        p = c;
        c1 = 2*p + 1;
        c2 = c1 + 1;
    }
}

int meshCheck(struct Mesh *m)
{
    size_t p;
    for(p=0; p<=(m->N-2)/2; p++)
    {
        size_t c1 = 2*p+1;
        size_t c2 = c1 + 1;

        if(c1 < m->N && m->heap[c1].err > m->heap[p].err)
            return 0;
        
        if(c2 < m->N && m->heap[c2].err > m->heap[p].err)
            return 0;

    }

    return 1;
}

void meshWrite(struct Mesh *m, char **buf)
{
    *buf = (char *) malloc((m->N * 4 * 30 + 12) * sizeof(char));

    size_t i;
    int c = sprintf(*buf, "%lu", m->N);
    for(i=0; i < m->N; i++)
    {
        struct Interval *in = &(m->heap[i]);
        c += sprintf(*buf + c, " %.16e %.16e %.16e %.16e",
                     in->a, in->b, in->I, in->err);
    }
    *buf = (char *) realloc(*buf, (c+1) * sizeof(char));
}

/***  Mesh 3  ***/

void mesh3Init(struct Mesh3 *m)
{
    size_t size = 4;
    m->totalSize = size;
    m->N = 0;
    m->heap = (struct Interval3 *)malloc(size * sizeof(struct Interval3));
}

void mesh3Free(struct Mesh3 *m)
{
    m->totalSize = 0;
    m->N = 0;
    free(m->heap);
    m->heap = NULL;
}

void mesh3Insert(struct Mesh3 *m, struct Interval3 *i)
{
    // Resize if necessary
    while(m->N >= m->totalSize)
    {
        m->totalSize *= 2;
        m->heap = (struct Interval3 *)realloc(m->heap,
                                    m->totalSize * sizeof(struct Interval3));
    }
    // Add interval to end of heap
    m->heap[m->N] = *i;
    (m->N)++;

    // Restore ordering
    mesh3HeapifyUp(m);
}

void mesh3Extract(struct Mesh3 *m, struct Interval3 *worst)
{
    *worst = m->heap[0];

    m->heap[0] = m->heap[m->N-1];
    (m->N)--;

    mesh3HeapifyDown(m);
}

double mesh3TotalIntegral(struct Mesh3 *m)
{
    double I = 0.0;

    size_t i = 0;
    for(i=0; i<m->N; i++)
        I += m->heap[i].I;
    
    return I;
}

double mesh3TotalError(struct Mesh3 *m)
{
    double err = 0.0;

    size_t i = 0;
    for(i=0; i<m->N; i++)
        err += m->heap[i].err;
    
    return err;
}

void mesh3HeapifyUp(struct Mesh3 *m)
{
    size_t c = m->N-1;
    size_t p = (c-1)/2;

    while(c != 0 && m->heap[p].err < m->heap[c].err)
    {
        struct Interval3 tempP = m->heap[p];
        m->heap[p] = m->heap[c];
        m->heap[c] = tempP;
        c = p;
        p = (c-1)/2;
    }
}

void mesh3HeapifyDown(struct Mesh3 *m)
{
    size_t p = 0;
    size_t c1 = 2*p + 1;
    size_t c2 = c1 + 1;


    while(c1 < m->N)
    {
        //Find the child with largest error
        size_t c = c1;
        double e = m->heap[c1].err;

        if(c2 < m->N && m->heap[c2].err > e)
        {
            c = c2;
            e = m->heap[c2].err;
        }

        // If the child is already in order then we're done.
        if(e <= m->heap[p].err)
            break;

        // Otherwise, swap with the child.
        struct Interval3 tempP = m->heap[p];
        m->heap[p] = m->heap[c];
        m->heap[c] = tempP;
        p = c;
        c1 = 2*p + 1;
        c2 = c1 + 1;
    }
}

int mesh3Check(struct Mesh3 *m)
{
    size_t p;
    for(p=0; p<=(m->N-2)/2; p++)
    {
        size_t c1 = 2*p+1;
        size_t c2 = c1 + 1;

        if(c1 < m->N && m->heap[c1].err > m->heap[p].err)
            return 0;
        
        if(c2 < m->N && m->heap[c2].err > m->heap[p].err)
            return 0;

    }

    return 1;
}

void mesh3Write(struct Mesh3 *m, char **buf)
{
    *buf = (char *) malloc((m->N * 4 * 30 + 12) * sizeof(char));

    size_t i;
    int c = sprintf(*buf, "%lu", m->N);
    for(i=0; i < m->N; i++)
    {
        struct Interval3 *in = &(m->heap[i]);
        c += sprintf(*buf + c, " %.16e %.16e %.16e %.16e",
                     in->a, in->b, in->I, in->err);
    }
    *buf = (char *) realloc(*buf, (c+1) * sizeof(char));
}

/***  Mesh 5  ***/

void mesh5Init(struct Mesh5 *m)
{
    size_t size = 4;
    m->totalSize = size;
    m->N = 0;
    m->heap = (struct Interval5 *)malloc(size * sizeof(struct Interval5));
}

void mesh5Free(struct Mesh5 *m)
{
    m->totalSize = 0;
    m->N = 0;
    free(m->heap);
    m->heap = NULL;
}

void mesh5Insert(struct Mesh5 *m, struct Interval5 *i)
{
    // Resize if necessary
    while(m->N >= m->totalSize)
    {
        m->totalSize *= 2;
        m->heap = (struct Interval5 *)realloc(m->heap,
                                    m->totalSize * sizeof(struct Interval5));
    }
    // Add interval to end of heap
    m->heap[m->N] = *i;
    (m->N)++;

    // Restore ordering
    mesh5HeapifyUp(m);
}

void mesh5Extract(struct Mesh5 *m, struct Interval5 *worst)
{
    *worst = m->heap[0];

    m->heap[0] = m->heap[m->N-1];
    (m->N)--;

    mesh5HeapifyDown(m);
}

double mesh5TotalIntegral(struct Mesh5 *m)
{
    double I = 0.0;

    size_t i = 0;
    for(i=0; i<m->N; i++)
        I += m->heap[i].I;
    
    return I;
}

double mesh5TotalError(struct Mesh5 *m)
{
    double err = 0.0;

    size_t i = 0;
    for(i=0; i<m->N; i++)
        err += m->heap[i].err;
    
    return err;
}

void mesh5HeapifyUp(struct Mesh5 *m)
{
    size_t c = m->N-1;
    size_t p = (c-1)/2;

    while(c != 0 && m->heap[p].err < m->heap[c].err)
    {
        struct Interval5 tempP = m->heap[p];
        m->heap[p] = m->heap[c];
        m->heap[c] = tempP;
        c = p;
        p = (c-1)/2;
    }
}

void mesh5HeapifyDown(struct Mesh5 *m)
{
    size_t p = 0;
    size_t c1 = 2*p + 1;
    size_t c2 = c1 + 1;


    while(c1 < m->N)
    {
        //Find the child with largest error
        size_t c = c1;
        double e = m->heap[c1].err;

        if(c2 < m->N && m->heap[c2].err > e)
        {
            c = c2;
            e = m->heap[c2].err;
        }

        // If the child is already in order then we're done.
        if(e <= m->heap[p].err)
            break;

        // Otherwise, swap with the child.
        struct Interval5 tempP = m->heap[p];
        m->heap[p] = m->heap[c];
        m->heap[c] = tempP;
        p = c;
        c1 = 2*p + 1;
        c2 = c1 + 1;
    }
}

int mesh5Check(struct Mesh5 *m)
{
    size_t p;
    for(p=0; p<=(m->N-2)/2; p++)
    {
        size_t c1 = 2*p+1;
        size_t c2 = c1 + 1;

        if(c1 < m->N && m->heap[c1].err > m->heap[p].err)
            return 0;
        
        if(c2 < m->N && m->heap[c2].err > m->heap[p].err)
            return 0;

    }

    return 1;
}

void mesh5Write(struct Mesh5 *m, char **buf)
{
    *buf = (char *) malloc((m->N * 4 * 30 + 12) * sizeof(char));

    size_t i;
    int c = sprintf(*buf, "%lu", m->N);
    for(i=0; i < m->N; i++)
    {
        struct Interval5 *in = &(m->heap[i]);
        c += sprintf(*buf + c, " %.16e %.16e %.16e %.16e",
                     in->a, in->b, in->I, in->err);
    }
    *buf = (char *) realloc(*buf, (c+1) * sizeof(char));
}

/***  Mesh 9  ***/

void mesh9Init(struct Mesh9 *m)
{
    size_t size = 4;
    m->totalSize = size;
    m->N = 0;
    m->heap = (struct Interval9 *)malloc(size * sizeof(struct Interval9));
}

void mesh9Free(struct Mesh9 *m)
{
    m->totalSize = 0;
    m->N = 0;
    free(m->heap);
    m->heap = NULL;
}

void mesh9Insert(struct Mesh9 *m, struct Interval9 *i)
{
    // Resize if necessary
    while(m->N >= m->totalSize)
    {
        m->totalSize *= 2;
        m->heap = (struct Interval9 *)realloc(m->heap,
                                    m->totalSize * sizeof(struct Interval9));
    }
    // Add interval to end of heap
    m->heap[m->N] = *i;
    (m->N)++;

    // Restore ordering
    mesh9HeapifyUp(m);
}

void mesh9Extract(struct Mesh9 *m, struct Interval9 *worst)
{
    *worst = m->heap[0];

    m->heap[0] = m->heap[m->N-1];
    (m->N)--;

    mesh9HeapifyDown(m);
}

double mesh9TotalIntegral(struct Mesh9 *m)
{
    double I = 0.0;

    size_t i = 0;
    for(i=0; i<m->N; i++)
        I += m->heap[i].I;
    
    return I;
}

double mesh9TotalError(struct Mesh9 *m)
{
    double err = 0.0;

    size_t i = 0;
    for(i=0; i<m->N; i++)
        err += m->heap[i].err;
    
    return err;
}

void mesh9HeapifyUp(struct Mesh9 *m)
{
    size_t c = m->N-1;
    size_t p = (c-1)/2;

    while(c != 0 && m->heap[p].err < m->heap[c].err)
    {
        struct Interval9 tempP = m->heap[p];
        m->heap[p] = m->heap[c];
        m->heap[c] = tempP;
        c = p;
        p = (c-1)/2;
    }
}

void mesh9HeapifyDown(struct Mesh9 *m)
{
    size_t p = 0;
    size_t c1 = 2*p + 1;
    size_t c2 = c1 + 1;


    while(c1 < m->N)
    {
        //Find the child with largest error
        size_t c = c1;
        double e = m->heap[c1].err;

        if(c2 < m->N && m->heap[c2].err > e)
        {
            c = c2;
            e = m->heap[c2].err;
        }

        // If the child is already in order then we're done.
        if(e <= m->heap[p].err)
            break;

        // Otherwise, swap with the child.
        struct Interval9 tempP = m->heap[p];
        m->heap[p] = m->heap[c];
        m->heap[c] = tempP;
        p = c;
        c1 = 2*p + 1;
        c2 = c1 + 1;
    }
}

int mesh9Check(struct Mesh9 *m)
{
    size_t p;
    if(m->N <= 1)
        return 1;

    for(p=0; p<=(m->N-2)/2; p++)
    {
        size_t c1 = 2*p+1;
        size_t c2 = c1 + 1;

        if(c1 < m->N && m->heap[c1].err > m->heap[p].err)
            return 0;
        
        if(c2 < m->N && m->heap[c2].err > m->heap[p].err)
            return 0;

    }

    return 1;
}

void mesh9Write(struct Mesh9 *m, char **buf)
{
    *buf = (char *) malloc((m->N * 4 * 30 + 12) * sizeof(char));

    size_t i;
    int c = sprintf(*buf, "%lu", m->N);
    for(i=0; i < m->N; i++)
    {
        struct Interval9 *in = &(m->heap[i]);
        c += sprintf(*buf + c, " %.16e %.16e %.16e %.16e",
                     in->a, in->b, in->I, in->err);
    }
    *buf = (char *) realloc(*buf, (c+1) * sizeof(char));
}

void interval9Write(struct Interval9 *i, FILE *stream)
{
    fprintf(stream, "(%.3le, %.3le)  %.12le +/- %.3le   %d\n",
            i->a, i->b, i->I, i->err, i->refinement);
    fprintf(stream, "   [%.3le %.3le %.3le %.3le %.3le %.3le"
                    " %.3le %.3le %.3le]\n", i->fa, i->fll, i->fl, i->flr,
                    i->fm, i->frl, i->fr, i->frr, i->fb);
}
