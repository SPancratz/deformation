/******************************************************************************

    Copyright (C) 2010, 2011 Sebastian Pancratz

******************************************************************************/

#ifndef STACK_H
#define STACK_H

#include <stdlib.h>

#define STACK_PROTOTYPE(NAME, TYPE, ATTR)                                     \
                                                                              \
typedef struct                                                                \
{                                                                             \
    TYPE * a;                                                                 \
    long n;                                                                   \
    long N;                                                                   \
} NAME ## _stack_struct;                                                      \
                                                                              \
typedef NAME ## _stack_struct NAME ## _stack_t[1];                            \
                                                                              \
ATTR void NAME ## _stack_init(NAME ## _stack_t S)                             \
{                                                                             \
    S->a = NULL;                                                              \
    S->n = 0;                                                                 \
    S->N = 0;                                                                 \
}                                                                             \
                                                                              \
ATTR void NAME ## _stack_init2(NAME ## _stack_t S, long N)                    \
{                                                                             \
    S->a = (TYPE *) malloc(N * sizeof(TYPE));                                 \
    S->n = 0;                                                                 \
    S->N = N;                                                                 \
}                                                                             \
                                                                              \
ATTR void NAME ## _stack_fit_size(NAME ## _stack_t S, long N)                 \
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
ATTR void NAME ## _stack_clear(NAME ## _stack_t S)                            \
{                                                                             \
    free(S->a);                                                               \
}                                                                             \
                                                                              \
ATTR int NAME ## _stack_is_empty(const NAME ## _stack_t S)                    \
{                                                                             \
    return S->n == 0;                                                         \
}                                                                             \
                                                                              \
ATTR void NAME ## _stack_push(NAME ## _stack_t S, TYPE o)                     \
{                                                                             \
    if (S->N == S->n)                                                         \
    {                                                                         \
        NAME ## _stack_fit_size(S, S->N + 1);                                 \
    }                                                                         \
    (S->a)[(S->n)++] = o;                                                     \
}                                                                             \
                                                                              \
ATTR TYPE NAME ## _stack_pop(NAME ## _stack_t S)                              \
{                                                                             \
    return (S->a)[--(S->n)];                                                  \
}                                                                             \

#endif

