/******************************************************************************

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#ifndef QUEUE_H
#define QUEUE_H

#include <stdlib.h>

#define QUEUE_PROTOTYPE(NAME, TYPE, ATTR)                                     \
                                                                              \
typedef struct                                                                \
{                                                                             \
    TYPE * a;                                                                 \
    long h;                                                                   \
    long n;                                                                   \
    long N;                                                                   \
} NAME ## _queue_struct;                                                      \
                                                                              \
typedef NAME ## _queue_struct NAME ## _queue_t[1];                            \
                                                                              \
ATTR void NAME ## _queue_init(NAME ## _queue_t Q)                             \
{                                                                             \
    Q->a = NULL;                                                              \
    Q->h = 0;                                                                 \
    Q->n = 0;                                                                 \
    Q->N = 0;                                                                 \
}                                                                             \
                                                                              \
ATTR void NAME ## _queue_init2(NAME ## _queue_t S, long N)                    \
{                                                                             \
    S->a = (TYPE *) malloc(N * sizeof(TYPE));                                 \
    S->h = 0;                                                                 \
    S->n = 0;                                                                 \
    S->N = N;                                                                 \
}                                                                             \
                                                                              \
ATTR void NAME ## _queue_fit_size(NAME ## _queue_t S, long N)                 \
{                                                                             \
    if (S->N == 0)                                                            \
    {                                                                         \
        S->a = (TYPE *) malloc(N * sizeof(TYPE));                             \
        S->N = N;                                                             \
    }                                                                         \
    else if (N > S->N)                                                        \
    {                                                                         \
        if (N < 2 * S->N)                                                     \
            N = 2 * S->N;                                                     \
        S->a = (TYPE *) realloc(S->a, N * sizeof(TYPE));                      \
        S->N = N;                                                             \
    }                                                                         \
}                                                                             \
                                                                              \
ATTR void NAME ## _queue_clear(NAME ## _queue_t S)                            \
{                                                                             \
    free(S->a);                                                               \
}                                                                             \
                                                                              \
ATTR int NAME ## _queue_is_empty(NAME ## _queue_t S)                          \
{                                                                             \
    return S->n == 0;                                                         \
}                                                                             \
                                                                              \
ATTR void NAME ## _queue_enqueue(NAME ## _queue_t S, TYPE o)                  \
{                                                                             \
    if (S->N == S->n)                                                         \
    {                                                                         \
        NAME ## _queue_fit_size(S, S->N + 1);                                 \
    }                                                                         \
    (S->a)[(S->h + (S->n)++) % S->N] = o;                                     \
}                                                                             \
                                                                              \
ATTR TYPE NAME ## _queue_dequeue(NAME ## _queue_t S)                          \
{                                                                             \
    long h = S->h;                                                            \
    (S->h) = (S->h + 1) % S->N;                                               \
    (S->n)--;                                                                 \
    return (S->a)[h];                                                         \
}                                                                             \

#endif

