/******************************************************************************

    Copyright (C) 2010, 2011 Sebastian Pancratz

******************************************************************************/

#ifndef RBTREE_H
#define RBTREE_H

#include <stdlib.h>

/* Macros independent of the parameters **************************************/

#define RBTREE_RED 0
#define RBTREE_BLACK 1

#define RBTREE_ROOT(t)  (*(t))

#define RBTREE_GRANDPARENT(n)  ((n)->parent->parent)

#define RBTREE_SIBLING(n) \
    (((n) == (n)->parent->left) ? (n)->parent->right : (n)->parent->left)

#define RBTREE_UNCLE(n)  (RBTREE_SIBLING((n)->parent))

#define RBTREE_COLOR(n)  ((n) ? (n)->color : RBTREE_BLACK)


/* Macros to be called instead of functions **********************************/

#define RBTREE_T(NAME)  NAME ## _rbtree_t

#define RBTREE_ITER_T(NAME)  NAME ## _rbtree_iter_t

#define RBTREE_ITER_INIT(NAME, iter, dict)                                    \
    NAME ## _rbtree_iter_init(iter, dict)

#define RBTREE_ITER_NEXT(NAME, iter)                                          \
    NAME ## _rbtree_iter_next(iter)

#define RBTREE_ITER_CLEAR(NAME, iter)                                         \
    NAME ## _rbtree_iter_clear(iter)

#define RBTREE_INIT(NAME, t)                                                  \
    NAME ## _rbtree_init(t)

#define RBTREE_CLEAR(NAME, t)                                                 \
    NAME ## _rbtree_clear(t)

#define RBTREE_SWAP(NAME, t1, t2)                                             \
    NAME ## _rbtree_swap(t1, t2)

#define RBTREE_SIZE(NAME, t)                                                  \
    NAME ## _rbtree_size(t)

#define RBTREE_FIND(NAME, k2, v2, t, k)                                       \
    NAME ## _rbtree_find(k2, v2, t, k)

#define RBTREE_INSERT(NAME, k2, v2, t, k, v)                                  \
    NAME ## _rbtree_insert(k2, v2, t, k, v)

#define RBTREE_DELETE(NAME, k2, v2, t, k)                                     \
    NAME ## _rbtree_delete(k2, v2, t, k)

#define RBTREE_IS_EMPTY(NAME, t)                                              \
    NAME ## _rbtree_is_empty(t)

#define RBTREE_PROTOTYPE_H(NAME, KTYPE, VTYPE, CMP, KCLEAR, VCLEAR, ATTR)     \
                                                                              \
typedef struct NAME ## _rbtree_node                                           \
{                                                                             \
    KTYPE key;                                                                \
    VTYPE val;                                                                \
    struct NAME ## _rbtree_node * left;                                       \
    struct NAME ## _rbtree_node * right;                                      \
    struct NAME ## _rbtree_node * parent;                                     \
    int color;                                                                \
} NAME ## _rbtree_node;                                                       \
                                                                              \
typedef NAME ## _rbtree_node * NAME ## _rbtree_t[1];                          \
                                                                              \
ATTR void                                                                     \
NAME ## _rbtree_init(NAME ## _rbtree_t t);                                    \
                                                                              \
ATTR void                                                                     \
NAME ## _rbtree_clear(NAME ## _rbtree_t t);                                   \
                                                                              \
ATTR void                                                                     \
NAME ## _rbtree_swap(NAME ## _rbtree_t t1, NAME ## _rbtree_t t2);             \
                                                                              \
long NAME ## _rbtree_size(const NAME ## _rbtree_t t);                         \
                                                                              \
ATTR NAME ## _rbtree_node *                                                   \
NAME ## _rbtree_find_node(const NAME ## _rbtree_t t, const KTYPE key);        \
                                                                              \
ATTR void                                                                     \
NAME ## _rbtree_replace_node(NAME ## _rbtree_t t,                             \
                             NAME ## _rbtree_node * o,                        \
                             NAME ## _rbtree_node * n);                       \
                                                                              \
ATTR void                                                                     \
NAME ## _rbtree_rotate_left(NAME ## _rbtree_t t,                              \
                            NAME ## _rbtree_node * n);                        \
                                                                              \
ATTR void                                                                     \
NAME ## _rbtree_rotate_right(NAME ## _rbtree_t t,                             \
                             NAME ## _rbtree_node * n);                       \
                                                                              \
ATTR int                                                                      \
NAME ## _rbtree_is_empty(const NAME ## _rbtree_t t);                          \
                                                                              \
ATTR NAME ## _rbtree_node *                                                   \
NAME ## _rbtree_min(const NAME ## _rbtree_node * n);                          \
                                                                              \
ATTR NAME ## _rbtree_node *                                                   \
NAME ## _rbtree_max(const NAME ## _rbtree_node * n);                          \
                                                                              \
ATTR NAME ## _rbtree_node *                                                   \
NAME ## _rbtree_prev(const NAME ## _rbtree_node * n);                         \
                                                                              \
ATTR NAME ## _rbtree_node *                                                   \
NAME ## _rbtree_next(const NAME ## _rbtree_node * n);                         \
                                                                              \
ATTR int                                                                      \
NAME ## _rbtree_find(KTYPE * key, VTYPE * val,                                \
                     const NAME ## _rbtree_t t, const KTYPE k);               \
                                                                              \
ATTR int                                                                      \
NAME ## _rbtree_insert(KTYPE * okey, VTYPE * oval,                            \
                       NAME ## _rbtree_t t, const KTYPE key, const VTYPE val);\
                                                                              \
ATTR int                                                                      \
NAME ## _rbtree_delete(KTYPE * okey, VTYPE * oval,                            \
                       NAME ## _rbtree_t t, const KTYPE key);                 \
                                                                              \
typedef struct NAME ## _rbtree_iter_struct                                    \
{                                                                             \
    NAME ## _rbtree_node ** S;                                                \
    long n;                                                                   \
} NAME ## _rbtree_iter_struct;                                                \
                                                                              \
typedef NAME ## _rbtree_iter_struct NAME ## _rbtree_iter_t[1];                \
                                                                              \
ATTR void                                                                     \
NAME ## _rbtree_iter_init(NAME ## _rbtree_iter_t iter,                        \
                          const NAME ## _rbtree_t t);                         \
                                                                              \
ATTR void                                                                     \
NAME ## _rbtree_iter_clear(NAME ## _rbtree_iter_t iter);                      \
                                                                              \
ATTR NAME ## _rbtree_node *                                                   \
NAME ## _rbtree_iter_next(NAME ## _rbtree_iter_t iter);                       \


#define RBTREE_PROTOTYPE_C(NAME, KTYPE, VTYPE, CMP, KCLEAR, VCLEAR, ATTR)     \
                                                                              \
ATTR void                                                                     \
NAME ## _rbtree_init(NAME ## _rbtree_t t)                                     \
{                                                                             \
    RBTREE_ROOT(t) = NULL;                                                    \
}                                                                             \
                                                                              \
ATTR void                                                                     \
_ ## NAME ## _rbtree_clear(NAME ## _rbtree_node * n)                          \
{                                                                             \
    if (n->left)                                                              \
        _ ## NAME ## _rbtree_clear(n->left);                                  \
    if (n->right)                                                             \
        _ ## NAME ## _rbtree_clear(n->right);                                 \
    KCLEAR(n->key);                                                           \
    VCLEAR(n->val);                                                           \
    free(n);                                                                  \
}                                                                             \
                                                                              \
ATTR void                                                                     \
NAME ## _rbtree_clear(NAME ## _rbtree_t t)                                    \
{                                                                             \
    if (RBTREE_ROOT(t))                                                       \
        _ ## NAME ## _rbtree_clear(RBTREE_ROOT(t));                           \
    RBTREE_ROOT(t) = NULL;                                                    \
}                                                                             \
                                                                              \
ATTR void                                                                     \
NAME ## _rbtree_swap(NAME ## _rbtree_t t1, NAME ## _rbtree_t t2)              \
{                                                                             \
    NAME ## _rbtree_node * n;                                                 \
    n = RBTREE_ROOT(t1);                                                      \
    RBTREE_ROOT(t1) = RBTREE_ROOT(t2);                                        \
    RBTREE_ROOT(t2) = n;                                                      \
}                                                                             \
                                                                              \
long NAME ## _rbtree_size(const NAME ## _rbtree_t t)                          \
{                                                                             \
    long N = 0;                                                               \
    NAME ## _rbtree_node * n;                                                 \
    NAME ## _rbtree_iter_t iter;                                              \
    NAME ## _rbtree_iter_init(iter, t);                                       \
    while ((n = NAME ## _rbtree_iter_next(iter)))                             \
        N++;                                                                  \
    NAME ## _rbtree_iter_clear(iter);                                         \
    return N;                                                                 \
}                                                                             \
                                                                              \
ATTR NAME ## _rbtree_node *                                                   \
NAME ## _rbtree_find_node(const NAME ## _rbtree_t t, const KTYPE key)         \
{                                                                             \
    NAME ## _rbtree_node * n = RBTREE_ROOT(t);                                \
    while (n)                                                                 \
    {                                                                         \
        int cmp = CMP(key, n->key);                                           \
        if (cmp == 0)                                                         \
            return n;                                                         \
        n = (cmp < 0) ? n->left : n->right;                                   \
    }                                                                         \
    return NULL;                                                              \
}                                                                             \
                                                                              \
ATTR void                                                                     \
NAME ## _rbtree_replace_node(NAME ## _rbtree_t t,                             \
                             NAME ## _rbtree_node * o,                        \
                             NAME ## _rbtree_node * n)                        \
{                                                                             \
    if (o->parent == NULL)                                                    \
    {                                                                         \
        RBTREE_ROOT(t) = n;                                                   \
    }                                                                         \
    else                                                                      \
    {                                                                         \
        if (o == o->parent->left)                                             \
            o->parent->left = n;                                              \
        else                                                                  \
            o->parent->right = n;                                             \
    }                                                                         \
    if (n)                                                                    \
    {                                                                         \
        n->parent = o->parent;                                                \
    }                                                                         \
}                                                                             \
                                                                              \
ATTR void                                                                     \
NAME ## _rbtree_rotate_left(NAME ## _rbtree_t t,                              \
                            NAME ## _rbtree_node * n)                         \
{                                                                             \
    NAME ## _rbtree_node * r = n->right;                                      \
                                                                              \
    NAME ## _rbtree_replace_node(t, n, r);                                    \
    n->right = r->left;                                                       \
    if (r->left)                                                              \
        r->left->parent = n;                                                  \
    r->left = n;                                                              \
    n->parent = r;                                                            \
}                                                                             \
                                                                              \
ATTR void                                                                     \
NAME ## _rbtree_rotate_right(NAME ## _rbtree_t t,                             \
                             NAME ## _rbtree_node * n)                        \
{                                                                             \
    NAME ## _rbtree_node *l = n->left;                                        \
                                                                              \
    NAME ## _rbtree_replace_node(t, n, l);                                    \
    n->left = l->right;                                                       \
    if (l->right)                                                             \
        l->right->parent = n;                                                 \
    l->right = n;                                                             \
    n->parent = l;                                                            \
}                                                                             \
                                                                              \
ATTR void                                                                     \
NAME ## _rbtree_insert_rearrange(NAME ## _rbtree_t t,                         \
                                 NAME ## _rbtree_node * n)                    \
{                                                                             \
    if (n->parent == NULL)  /* Case 1 */                                      \
        n->color = RBTREE_BLACK;                                              \
                                                                              \
    else if (n->parent->color == RBTREE_BLACK)  /* Case 2 */                  \
        return;                                                               \
                                                                              \
    else if (RBTREE_COLOR(RBTREE_UNCLE(n)) == RBTREE_RED)  /* Case 3 */       \
    {                                                                         \
        n->parent->color = RBTREE_BLACK;                                      \
        RBTREE_UNCLE(n)->color = RBTREE_BLACK;                                \
        RBTREE_GRANDPARENT(n)->color = RBTREE_RED;                            \
        NAME ## _rbtree_insert_rearrange(t, RBTREE_GRANDPARENT(n));           \
    }                                                                         \
    else  /* Cases 4, 5 */                                                    \
    {                                                                         \
        if (n == n->parent->right && n->parent == RBTREE_GRANDPARENT(n)->left)\
        {                                                                     \
            NAME ## _rbtree_rotate_left(t, n->parent);                        \
            n = n->left;                                                      \
        }                                                                     \
        else if (n == n->parent->left && n->parent == RBTREE_GRANDPARENT(n)->right) \
        {                                                                     \
            NAME ## _rbtree_rotate_right(t, n->parent);                       \
            n = n->right;                                                     \
        }                                                                     \
                                                                              \
        n->parent->color = RBTREE_BLACK;                                      \
        RBTREE_GRANDPARENT(n)->color = RBTREE_RED;                            \
        if (n == n->parent->left && n->parent == RBTREE_GRANDPARENT(n)->left) \
            NAME ## _rbtree_rotate_right(t, RBTREE_GRANDPARENT(n));           \
        else                                                                  \
            NAME ## _rbtree_rotate_left(t, RBTREE_GRANDPARENT(n));            \
    }                                                                         \
}                                                                             \
                                                                              \
ATTR void                                                                     \
NAME ## _rbtree_delete_rearrange(NAME ## _rbtree_t t,                         \
                                 NAME ## _rbtree_node * n)                    \
{                                                                             \
    if (n->parent == NULL)                                                    \
        return;                                                               \
                                                                              \
    if (RBTREE_COLOR(RBTREE_SIBLING(n)) == RBTREE_RED)                        \
    {                                                                         \
        n->parent->color = RBTREE_RED;                                        \
        RBTREE_SIBLING(n)->color = RBTREE_BLACK;                              \
        if (n == n->parent->left)                                             \
            NAME ## _rbtree_rotate_left(t, n->parent);                        \
        else                                                                  \
            NAME ## _rbtree_rotate_right(t, n->parent);                       \
    }                                                                         \
                                                                              \
    if (RBTREE_COLOR(n->parent) == RBTREE_BLACK &&                            \
        RBTREE_COLOR(RBTREE_SIBLING(n)) == RBTREE_BLACK &&                    \
        RBTREE_COLOR(RBTREE_SIBLING(n)->left) == RBTREE_BLACK &&              \
        RBTREE_COLOR(RBTREE_SIBLING(n)->right) == RBTREE_BLACK)               \
    {                                                                         \
        RBTREE_SIBLING(n)->color = RBTREE_RED;                                \
        NAME ## _rbtree_delete_rearrange(t, n->parent);                       \
    }                                                                         \
    else if (RBTREE_COLOR(n->parent) == RBTREE_RED &&                         \
             RBTREE_COLOR(RBTREE_SIBLING(n)) == RBTREE_BLACK &&               \
             RBTREE_COLOR(RBTREE_SIBLING(n)->left) == RBTREE_BLACK &&         \
             RBTREE_COLOR(RBTREE_SIBLING(n)->right) == RBTREE_BLACK)          \
    {                                                                         \
        RBTREE_SIBLING(n)->color = RBTREE_RED;                                \
        n->parent->color = RBTREE_BLACK;                                      \
    }                                                                         \
    else                                                                      \
    {                                                                         \
        if (n == n->parent->left &&                                           \
            RBTREE_COLOR(RBTREE_SIBLING(n)) == RBTREE_BLACK &&                \
            RBTREE_COLOR(RBTREE_SIBLING(n)->left) == RBTREE_RED &&            \
            RBTREE_COLOR(RBTREE_SIBLING(n)->right) == RBTREE_BLACK)           \
        {                                                                     \
            RBTREE_SIBLING(n)->color = RBTREE_RED;                            \
            RBTREE_SIBLING(n)->left->color = RBTREE_BLACK;                    \
            NAME ## _rbtree_rotate_right(t, RBTREE_SIBLING(n));               \
        }                                                                     \
        else if (n == n->parent->right &&                                     \
                 RBTREE_COLOR(RBTREE_SIBLING(n)) == RBTREE_BLACK &&           \
                 RBTREE_COLOR(RBTREE_SIBLING(n)->right) == RBTREE_RED &&      \
                 RBTREE_COLOR(RBTREE_SIBLING(n)->left) == RBTREE_BLACK)       \
        {                                                                     \
            RBTREE_SIBLING(n)->color = RBTREE_RED;                            \
            RBTREE_SIBLING(n)->right->color = RBTREE_BLACK;                   \
            NAME ## _rbtree_rotate_left(t, RBTREE_SIBLING(n));                \
        }                                                                     \
                                                                              \
        RBTREE_SIBLING(n)->color = RBTREE_COLOR(n->parent);                   \
        n->parent->color = RBTREE_BLACK;                                      \
        if (n == n->parent->left)                                             \
        {                                                                     \
            RBTREE_SIBLING(n)->right->color = RBTREE_BLACK;                   \
            NAME ## _rbtree_rotate_left(t, n->parent);                        \
        }                                                                     \
        else                                                                  \
        {                                                                     \
            RBTREE_SIBLING(n)->left->color = RBTREE_BLACK;                    \
            NAME ## _rbtree_rotate_right(t, n->parent);                       \
        }                                                                     \
    }                                                                         \
}                                                                             \
                                                                              \
ATTR NAME ## _rbtree_node *                                                   \
NAME ## _rbtree_new_node(KTYPE key, VTYPE val,                                \
                         int color,                                           \
                         NAME ## _rbtree_node * parent,                       \
                         NAME ## _rbtree_node * left,                         \
                         NAME ## _rbtree_node * right)                        \
{                                                                             \
    NAME ## _rbtree_node * n;                                                 \
                                                                              \
    n = malloc(sizeof(struct NAME ## _rbtree_node));                          \
    n->key = key;                                                             \
    n->val = val;                                                             \
    n->color = color;                                                         \
    n->parent = parent;                                                       \
    n->left = left;                                                           \
    n->right = right;                                                         \
    if (left) n->left->parent = n;                                            \
    if (right) n->right->parent = n;                                          \
    return n;                                                                 \
}                                                                             \
                                                                              \
ATTR int                                                                      \
NAME ## _rbtree_is_empty(const NAME ## _rbtree_t t)                           \
{                                                                             \
    return RBTREE_ROOT(t) == NULL;                                            \
}                                                                             \
                                                                              \
ATTR NAME ## _rbtree_node *                                                   \
NAME ## _rbtree_min(const NAME ## _rbtree_node * n)                           \
{                                                                             \
    NAME ## _rbtree_node *m = (NAME ## _rbtree_node *) n;                     \
    while (m->left)                                                           \
        m = m->left;                                                          \
    return m;                                                                 \
}                                                                             \
                                                                              \
ATTR NAME ## _rbtree_node *                                                   \
NAME ## _rbtree_max(const NAME ## _rbtree_node * n)                           \
{                                                                             \
    NAME ## _rbtree_node *m = (NAME ## _rbtree_node *) n;                     \
    while (m->right)                                                          \
        m = m->right;                                                         \
    return m;                                                                 \
}                                                                             \
                                                                              \
ATTR NAME ## _rbtree_node *                                                   \
NAME ## _rbtree_prev(const NAME ## _rbtree_node * n)                          \
{                                                                             \
    if (n->left)                                                              \
    {                                                                         \
        return NAME ## _rbtree_max(n->left);                                  \
    }                                                                         \
    else                                                                      \
    {                                                                         \
        NAME ## _rbtree_node *y, *z;                                          \
                                                                              \
        y = (NAME ## _rbtree_node *) n;                                       \
        z = (NAME ## _rbtree_node *) n->parent;                               \
                                                                              \
        while ((z) && y == z->left)                                           \
        {                                                                     \
            y = z;                                                            \
            z = z->parent;                                                    \
        }                                                                     \
        return z;                                                             \
    }                                                                         \
}                                                                             \
                                                                              \
ATTR NAME ## _rbtree_node *                                                   \
NAME ## _rbtree_next(const NAME ## _rbtree_node * n)                          \
{                                                                             \
    if (n->right)                                                             \
    {                                                                         \
        return NAME ## _rbtree_min(n->right);                                 \
    }                                                                         \
    else                                                                      \
    {                                                                         \
        NAME ## _rbtree_node *y, *z;                                          \
                                                                              \
        y = (NAME ## _rbtree_node *) n;                                       \
        z = (NAME ## _rbtree_node *) n->parent;                               \
                                                                              \
        while ((z) && y == z->right)                                          \
        {                                                                     \
            y = z;                                                            \
            z = z->parent;                                                    \
        }                                                                     \
        return z;                                                             \
    }                                                                         \
}                                                                             \
                                                                              \
ATTR int                                                                      \
NAME ## _rbtree_find(KTYPE * key, VTYPE * val,                                \
                     const NAME ## _rbtree_t t, const KTYPE k)                \
{                                                                             \
    NAME ## _rbtree_node * n = NAME ## _rbtree_find_node(t, k);               \
                                                                              \
    if (n)                                                                    \
    {                                                                         \
        *key = n->key;                                                        \
        *val = n->val;                                                        \
        return 1;                                                             \
    }                                                                         \
    return 0;                                                                 \
}                                                                             \
                                                                              \
ATTR int                                                                      \
NAME ## _rbtree_insert(KTYPE * okey, VTYPE * oval,                            \
                       NAME ## _rbtree_t t, const KTYPE key, const VTYPE val) \
{                                                                             \
    if (RBTREE_ROOT(t) == NULL)                                               \
    {                                                                         \
        RBTREE_ROOT(t) =                                                      \
          NAME ## _rbtree_new_node(key, val, RBTREE_BLACK, NULL, NULL, NULL); \
        return 0;                                                             \
    }                                                                         \
    else                                                                      \
    {                                                                         \
        NAME ## _rbtree_node * n = RBTREE_ROOT(t);                            \
        int cmp = CMP(key, n->key);                                           \
                                                                              \
        while (cmp)                                                           \
        {                                                                     \
            if (cmp < 0)                                                      \
            {                                                                 \
                if (n->left)                                                  \
                    n = n->left;                                              \
                else                                                          \
                {                                                             \
                    n->left = NAME ## _rbtree_new_node(key, val, RBTREE_RED,  \
                                                       n, NULL, NULL);        \
                    NAME ## _rbtree_insert_rearrange(t, n->left);             \
                    return 0;                                                 \
                }                                                             \
            }                                                                 \
            else                                                              \
            {                                                                 \
                if (n->right)                                                 \
                    n = n->right;                                             \
                else                                                          \
                {                                                             \
                    n->right = NAME ## _rbtree_new_node(key, val, RBTREE_RED, \
                                                        n, NULL, NULL);       \
                    NAME ## _rbtree_insert_rearrange(t, n->right);            \
                    return 0;                                                 \
                }                                                             \
            }                                                                 \
            cmp = CMP(key, n->key);                                           \
        }                                                                     \
                                                                              \
        *okey = n->key;                                                       \
        *oval = n->val;                                                       \
        n->key = key;                                                         \
        n->val = val;                                                         \
        return 1;                                                             \
    }                                                                         \
}                                                                             \
                                                                              \
ATTR int                                                                      \
NAME ## _rbtree_delete(KTYPE * okey, VTYPE * oval,                            \
                       NAME ## _rbtree_t t, const KTYPE key)                  \
{                                                                             \
    NAME ## _rbtree_node *c, *n = NAME ## _rbtree_find_node(t, key);          \
                                                                              \
    if (n == NULL)                                                            \
        return 0;                                                             \
                                                                              \
    *okey = n->key;                                                           \
    *oval = n->val;                                                           \
                                                                              \
    if ((n->left) && (n->right))                                              \
    {                                                                         \
        NAME ## _rbtree_node *p = n->left;                                    \
        while (p->right)                                                      \
            p = p->right;                                                     \
        n->key = p->key;                                                      \
        n->val = p->val;                                                      \
        n = p;                                                                \
    }                                                                         \
                                                                              \
    c = (n->right) ? n->right : n->left;                                      \
    if (n->color == RBTREE_BLACK)                                             \
    {                                                                         \
        n->color = RBTREE_COLOR(c);                                           \
        NAME ## _rbtree_delete_rearrange(t, n);                               \
    }                                                                         \
    NAME ## _rbtree_replace_node(t, n, c);                                    \
    if (n->parent == NULL && (c))                                             \
        c->color = RBTREE_BLACK;                                              \
    free(n);                                                                  \
    return 1;                                                                 \
}                                                                             \
                                                                              \
ATTR void                                                                     \
NAME ## _rbtree_iter_init(NAME ## _rbtree_iter_t iter,                        \
                          const NAME ## _rbtree_t t)                          \
{                                                                             \
    long N;                                                                   \
    NAME ## _rbtree_node * n;                                                 \
                                                                              \
    N = 0;                                                                    \
    n = RBTREE_ROOT(t);                                                       \
    while (n)                                                                 \
    {                                                                         \
        if (n->color == RBTREE_BLACK)                                         \
            N++;                                                              \
        n = n->left;                                                          \
    }                                                                         \
    N = 2 * N;                                                                \
                                                                              \
    iter->S = malloc(N * sizeof(NAME ## _rbtree_node *));                     \
    iter->n = 0;                                                              \
                                                                              \
    n = RBTREE_ROOT(t);                                                       \
    while (n)                                                                 \
    {                                                                         \
        iter->S[(iter->n)++] = n;                                             \
        n = n->left;                                                          \
    }                                                                         \
}                                                                             \
                                                                              \
ATTR void                                                                     \
NAME ## _rbtree_iter_clear(NAME ## _rbtree_iter_t iter)                       \
{                                                                             \
    free(iter->S);                                                            \
}                                                                             \
                                                                              \
ATTR NAME ## _rbtree_node *                                                   \
NAME ## _rbtree_iter_next(NAME ## _rbtree_iter_t iter)                        \
{                                                                             \
    if (iter->n)                                                              \
    {                                                                         \
        NAME ## _rbtree_node *n, *r = iter->S[--(iter->n)];                   \
                                                                              \
        n = r->right;                                                         \
        while (n)                                                             \
        {                                                                     \
            iter->S[(iter->n)++] = n;                                         \
            n = n->left;                                                      \
        }                                                                     \
                                                                              \
        return r;                                                             \
    }                                                                         \
    return NULL;                                                              \
}                                                                             \

/**
 * \def    RBTREE_PROTOTYPE_DEBUG_C(NAME, ATTR)
 * \brief  Generates code useful for debugging.
 * 
 * This macro prototypes the functions
 * 
 *     -# <code>NAME_rbtree_verify2()</code>
 *     -# <code>NAME_rbtree_verify4()</code>
 *     -# <code>NAME_rbtree_verify5()</code>
 * 
 * In each case, these functions return <code>1</code> if the property 
 * in question has been verified and <code>0</code> otherwise.
 */

#define RBTREE_PROTOTYPE_DEBUG_H(NAME, ATTR)                                  \
                                                                              \
ATTR int                                                                      \
NAME ## _rbtree_verify2(const NAME ## _rbtree_node * n);                      \
                                                                              \
ATTR int                                                                      \
NAME ## _rbtree_verify4(const NAME ## _rbtree_node * n);                      \
                                                                              \
ATTR int                                                                      \
NAME ## _rbtree_verify5(const NAME ## _rbtree_node * n);                      \


#define RBTREE_PROTOTYPE_DEBUG_C(NAME, ATTR)                                  \
                                                                              \
ATTR int                                                                      \
NAME ## _rbtree_verify2(const NAME ## _rbtree_node * n)                       \
{                                                                             \
    return (n == NULL || n->color == RBTREE_BLACK);                           \
}                                                                             \
                                                                              \
ATTR int                                                                      \
NAME ## _rbtree_verify4(const NAME ## _rbtree_node * n)                       \
{                                                                             \
    int ans;                                                                  \
                                                                              \
    if (n == NULL)                                                            \
        return 1;                                                             \
                                                                              \
    ans = (n->color == RBTREE_RED) ? n->parent->color == RBTREE_BLACK : 1;    \
                                                                              \
    ans = ans && NAME ## _rbtree_verify4(n->left);                            \
    ans = ans && NAME ## _rbtree_verify4(n->right);                           \
    return ans;                                                               \
}                                                                             \
                                                                              \
ATTR int                                                                      \
_ ## NAME ## _rbtree_verify5(const NAME ## _rbtree_node * n,                  \
                             int bc, int * pbc)                               \
{                                                                             \
    if (RBTREE_COLOR(n) == RBTREE_BLACK)                                      \
        bc++;                                                                 \
                                                                              \
    if (n == NULL)                                                            \
    {                                                                         \
        if (*pbc == -1)                                                       \
        {                                                                     \
            *pbc = bc;                                                        \
            return 1;                                                         \
        }                                                                     \
        return (bc == *pbc);                                                  \
    }                                                                         \
    else                                                                      \
    {                                                                         \
        int ans = 1;                                                          \
        ans = ans && _ ## NAME ## _rbtree_verify5(n->left, bc, pbc);          \
        ans = ans && _ ## NAME ## _rbtree_verify5(n->right, bc, pbc);         \
        return ans;                                                           \
    }                                                                         \
}                                                                             \
                                                                              \
ATTR int                                                                      \
NAME ## _rbtree_verify5(const NAME ## _rbtree_node * n)                       \
{                                                                             \
    int pbc = -1;                                                             \
    return _ ## NAME ## _rbtree_verify5(n, 0, &pbc);                          \
}                                                                             \

#endif

