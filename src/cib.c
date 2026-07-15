/*
 * CIB: Cofinite Integral Braces in GAP
 */

#include <gap_all.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#if defined(GAP_KERNEL_DEBUG) && GAP_KERNEL_DEBUG
#ifndef CIB_DIAG
#   define CIB_DIAG 1
#endif
#else
#ifndef CIB_DIAG
#   define CIB_DIAG 0
#endif
#endif

/* GAP type used in the shell */
static Obj CIBVectorSystemContextDataType;

/* Check for proper type of the context argument */
#define CHECK_ARGUMENT_CTX(ctx) \
    RequireArgumentCondition( SELF_NAME, ctx, IS_DATOBJ(ctx) && TYPE_OBJ(ctx) == CIBVectorSystemContextDataType, \
    "must be a CIB vector system context data object (T_DATOBJ)" )

/* aligning */
static unsigned VEC_BITS   = 128; // defaults to SSE2
static unsigned VEC_BYTES  = 16;  // defaults to SSE2
static unsigned HEAD_BYTES = 16;  // default space for beginning of the data bag
static unsigned HEAD_SLOTS = 16 / sizeof(Obj);

static inline unsigned detect_vec_bits_runtime(void) {
#if defined(__GNUC__) || defined(__clang__)
    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx512f")) return 512;
    if (__builtin_cpu_supports("avx2"))    return 256;
    return 128; // SSE2 fallback
#else
    return 128;
#endif
}

/* helpers */

/*
 * @brief Calculates the lowest multiple of n, greater than or equal to m
 *
 * @param m an integer
 * @param n an integer
 *
 * @return ⌈m/n⌉·n
 */
static inline size_t find_multiplicity_roof(const size_t m, const size_t n)
{
    GAP_ASSERT( n>0 );
    return ((m + n - 1) / n) * n;
}

/*
 * @brief next power of two
 *
 * Calculate the least power of two greater than or equal to x
 */
static inline size_t next_pow2(size_t x) {
    if (x <= 2) return 2;
    x--;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
#if SIZE_MAX > 0xFFFFFFFFu
    x |= x >> 32;
#endif
    if (x == SIZE_MAX) {
        /* mark overflow */
        return 0;
    }
    return x + 1;
}

/*
 * @brief Compute the dimension used to include SSE/AVX cpu instructions
 *
 * In trying to help the compiler to produce SSE/AVX friendly code,
 * we alter the dimension given by the user to the multiplicity of
 * bits used by the most efficient CPU instructions available.
 *
 */
static inline size_t compute_aligned_dimension(const size_t dim)
{
    const size_t lanes = VEC_BYTES / sizeof(uint16_t);
    if ( dim > lanes ) {
        return find_multiplicity_roof( dim, lanes );
    }
    return next_pow2(dim);
}

/* ========================== structure of the data bag ===============================
 *
 * Bag T_DATOBJ of size: HEAD_BYTES + total_bytes
 *
 * N = HEAD_SLOTS - 1
 *
 * [ slot0(T_DATOBJ) slot1(0) ... slotN(0) |      cib_data_t    | raw-matrices | tuple ]
 * [ <----------- HEAD_BYTES ------------> | <-- raw_offset --> |              |       ]
 * [                                       | <---------  tup_offset ---------> |       ]
 * [                                       | <-------------  total_bytes ------------> ]
 */

typedef struct {
    size_t   dim;          // dimension of the group
    size_t   dim_aligned;  // aligned dimension
    uint16_t exp;          // current exponent
    size_t   rank;         // rank of the coboundary group
    size_t   size;         // group order
    size_t   len;          // aligned length of one linearized matrix
    size_t   raw_offset;   // offset (od base) do danych macierzy
    size_t   tup_offset;   // offset (od base) do krotki współrzędnych
} cib_data_t;

/* raw pointer to the beginning of the data */
static inline uint8_t* INLINE_BASE(Obj ctx) {
    return (uint8_t*)ADDR_OBJ(ctx) + HEAD_BYTES;
}

/* size of the actual data */
static inline size_t INLINE_BYTES(Obj ctx) {
    size_t sz = SIZE_OBJ(ctx);
    return sz >= HEAD_BYTES ? (sz - HEAD_BYTES) : 0;
}

/* pointer to the data */
static inline cib_data_t *DATA_PTR_INLINE(Obj ctx) {
    return (cib_data_t*) INLINE_BASE(ctx);
}

/* pointer to the matrix part of the data */
static inline uint16_t *RAW_PTR(cib_data_t* data) {
    return (uint16_t*)((uint8_t*)data + data->raw_offset);
}
/* pointer to the specific matrix */
static inline uint16_t *IND_PTR(cib_data_t *data, const size_t ind) {
    return RAW_PTR(data) + ind * data->len;
}
/* pointer to the tuple */
static inline uint16_t *TUP_PTR(cib_data_t* data) {
    return (uint16_t*)((uint8_t*)data + data->tup_offset);
}

/* print some debug info */
static inline void dbg_dump_ctx(cib_data_t *d, const char *tag) {
#if defined(GAP_KERNEL_DEBUG) && GAP_KERNEL_DEBUG
    uint8_t *base = (uint8_t*)d;
    fprintf(stderr,
        "[cib:%s] ADDR=%08lx VEC_BYTES=%u size(rows)=%zu dim=%zu dim_aligned=%zu "
        "len=%zu rows*stride=%zu rank=%zu raw_offset=%lu tup_offset=%lu\n",
        tag, (uintptr_t)d, VEC_BYTES,
        d->size, d->dim, d->dim_aligned,
        d->len, d->size * d->dim_aligned, d->rank,
        d->raw_offset, d->tup_offset);
#else
    (void)d; (void)tag;
#endif
}

/*
 * @brief Convert GAP list to data in GAP bag
 *
 * @param list a GAP list
 * @param ctx  the GAP bag holding the data
 * @param ind  the number of the matrix to be filled in
 *
 * @return true if all elements of the list are integers, false otherwise
 */
static inline bool list_to_data(Obj list, Obj ctx, const size_t ind)
{
    cib_data_t    *data = DATA_PTR_INLINE(ctx);
    const size_t   size = data->size;
    const size_t   dim  = data->dim;
    const uint16_t exp  = data->exp;
    const size_t   pad  = data->dim_aligned - dim;
#if defined(GAP_KERNEL_DEBUG) && GAP_KERNEL_DEBUG
    const size_t   len  = data->len;
#endif
    uint16_t      *dst  = NULL; //IND_PTR(data, ind);

    size_t k = 0;
    for (size_t i=1; i<=size; ++i) {
        Obj row = ELM_PLIST(list, i);
        if (!IS_PLIST(row)) {
            ErrorQuit("list_to_data: list elements must be lists", 0, 0);
        }
        for (size_t j=1; j<=dim; j++) {
            Obj val = ELM_PLIST(row, j);
            if (!IS_INTOBJ(val)) {
                return false;
            }
            int num = INT_INTOBJ(val) % exp;
            /* possible gc */
            dst    = IND_PTR( DATA_PTR_INLINE(ctx), ind );
            dst[k] = (uint16_t)(num>=0 ? num : num + exp);
            k++;
            GAP_ASSERT( k <= len );
        }
        k += pad;
        GAP_ASSERT( k <= len );
    }
    return true;
}

/*
 * @brief Create GAP bag for calculations
 *
 * @param exp   the exponent - all calculations are made mod exp
 * @param vs    the vector system of the group
 * @param basis the basis of the 1-coboundary group of the holonomy/point group
 *
 * What is assumed is that vs and the elements of basis are lists with the property
 * that i-th elements on every one of them correspond to the same element of the point group.
 *
 * @return the GAP bag with the structure as described above.
 */
Obj FuncCIBVectorSystemContextCreate(Obj self, Obj exp, Obj vs, Obj basis)
{
    RequirePositiveSmallInt( SELF_NAME, exp );
    RequirePlainList( SELF_NAME, vs );
    RequirePlainList( SELF_NAME, basis );

    const size_t size = LEN_PLIST(vs);
    Obj row = ELM_PLIST(vs, 1);
    if (!IS_PLIST(row)) {
        ErrorQuit("CIBVectorSystemContextCreate: <vs> elements must be lists", 0, 0);
    }
    const size_t dim         = LEN_PLIST(row);
    const size_t dim_aligned = dim;             // in case we would like to use AVX aligning for rows of the matrices
    const size_t len         = dim * size;
    const size_t rank        = LEN_PLIST(basis);
    GAP_ASSERT( dim <= dim_aligned );

    /* calculate field sizes */
    const size_t cib_bytes   = find_multiplicity_roof(                    sizeof(cib_data_t), VEC_BYTES );
    const size_t raw_bytes   = find_multiplicity_roof( len * (rank + 1) * sizeof( uint16_t ), VEC_BYTES );
    const size_t tup_bytes   = find_multiplicity_roof(       (rank + 1) * sizeof( uint16_t ), VEC_BYTES );
    const size_t total_bytes = cib_bytes + raw_bytes + tup_bytes;

    Obj ctx = NewBag(T_DATOBJ, HEAD_BYTES + total_bytes);
    SET_TYPE_OBJ(ctx, CIBVectorSystemContextDataType);
    for (unsigned k=1; k<HEAD_SLOTS; ++k) ADDR_OBJ(ctx)[k] = (Obj)0;
    CHANGED_BAG(ctx);

    /* fill the data with zeroes */
    uint8_t *base = INLINE_BASE(ctx);
    memset(base, 0, total_bytes);

    cib_data_t *data = (cib_data_t*) base;

    data->exp         = (uint16_t)INT_INTOBJ(exp);
    data->size        = size;
    data->dim         = dim;
    data->dim_aligned = dim_aligned;
    data->rank        = rank;
    data->len         = len;
    data->raw_offset  = cib_bytes;
    data->tup_offset  = cib_bytes + raw_bytes;

    dbg_dump_ctx(data, "create1");

    if (!list_to_data( vs, ctx, 0 )) {
        return Fail;
    }
    CHANGED_BAG( ctx );
    for (size_t i=1; i<=rank; i++) {
        Obj vec = ELM_PLIST(basis, i);
        if (!IS_PLIST(vec)) {
            ErrorQuit("CIBVectorSystemContextCreate: basis elements must be lists", 0, 0);
        }
        if (LEN_PLIST(vec) != size) {
            ErrorQuit("CIBVectorSystemContextCreate: basis elements must have the same size as <vs>", 0, 0);
        }
        if (!list_to_data( vec, ctx, i )) {
            ErrorQuit("CIBVectorSystemContextCreate: invalid basis vector no %d", (Int)i, 0);
        }
        CHANGED_BAG( ctx );
    }

    dbg_dump_ctx( DATA_PTR_INLINE(ctx) , "create2" );
    return ctx;
}

static inline Obj obj_coordinates(Obj ctx)
{
    cib_data_t  *data = DATA_PTR_INLINE(ctx);
    const size_t rank = data->rank;

    Obj vec = NEW_PLIST(T_PLIST, rank);
    SET_LEN_PLIST(vec, rank);
    CHANGED_BAG( vec );

    for (size_t i=1; i<=rank; i++) {
        SET_ELM_PLIST(vec, i, ObjInt_UInt( TUP_PTR(DATA_PTR_INLINE(ctx))[i] ));
        CHANGED_BAG(vec);
    }
    return vec;
}

Obj FuncCIBVectorSystemContextRawCoordinates(Obj self, Obj ctx)
{
  CHECK_ARGUMENT_CTX(ctx);
  return obj_coordinates( ctx );
}

static inline void add_pos_mul(cib_data_t *data, const size_t ind, uint16_t val)
{
    GAP_ASSERT( ind>0 );

    uint16_t * __restrict src = IND_PTR(data, ind);
    uint16_t * __restrict dst = RAW_PTR(data);

    const uint64_t exp = data->exp;
    const size_t   len = data->len;

    for (size_t i=0; i<len; ++i) {
        uint64_t s = (uint64_t)dst[i] + (uint64_t)val * (uint64_t)src[i];
        dst[i] = (uint16_t)(s % exp);
    }
}

Obj FuncCIBVectorSystemContextRawSet(Obj self, Obj ctx, Obj coordinates)
{
    CHECK_ARGUMENT_CTX(ctx);

    RequirePossList( SELF_NAME, coordinates );

    cib_data_t *data = DATA_PTR_INLINE(ctx);
    uint16_t      *c = NULL; //TUP_PTR(data);
    size_t      rank = data->rank;
    uint16_t     exp = data->exp;

    // reset overflow marker
    c[0] = 0;
    // set the output vectors to zero
    memset( RAW_PTR(data), 0, data->len * sizeof(uint16_t) );

    for (size_t i=1; i<=rank; i++) {
        /* Does reading may trigger gc?  */
        Obj val  = ELM_PLIST(coordinates, i);

        int num  = INT_INTOBJ(val) % exp;
            data = DATA_PTR_INLINE(ctx);
            c    = TUP_PTR( data );
            c[i] = (uint16_t)(num >= 0 ? num : num + exp);

        add_pos_mul( data, i, c[i] );

        CHANGED_BAG( ctx );
    }
    return obj_coordinates( ctx );
}

static inline void add_pos(Obj ctx, const size_t ind)
{
    GAP_ASSERT( ind>0 );
    cib_data_t  *data  = DATA_PTR_INLINE( ctx );

#if defined(CIB_DIAG) && CIB_DIAG
    const size_t avail = INLINE_BYTES(ctx);
    const size_t need_raw = data->raw_offset + (data->rank + 1) * data->len * sizeof(uint16_t);
    if (need_raw > avail) {
        ErrorQuit("CIB/add_pos: raw region OOB (need=%d, have=%d)",
                (Int)need_raw, (Int)avail);
    }
    if (ind > data->rank) {
        ErrorQuit("CIB/add_pos: index %d out of range 1..%d",
                (Int)ind, (Int)data->rank);
    }
#endif

    uint16_t * __restrict src = IND_PTR(data, ind);
    uint16_t * __restrict dst = RAW_PTR(data);
    const uint32_t        exp = data->exp;
    const size_t          len = data->len;

    for (size_t i=0; i<len; ++i) {
        uint32_t s = (uint32_t)dst[i] + (uint32_t)src[i];
        dst[i] = (uint16_t)(s >= exp ? s - exp : s);
    }
    // not needed as we do not do anything to fire up gc?
    // CHANGED_BAG(ctx);
}

static inline bool inc(Obj ctx)
{
    cib_data_t *data = DATA_PTR_INLINE( ctx );

#if defined(CIB_DIAG) && CIB_DIAG
    const size_t avail = INLINE_BYTES(ctx);

    const size_t need_tup = data->tup_offset + (data->rank + 1) * sizeof(uint16_t);
    if (need_tup > avail) {
        ErrorQuit("CIB/inc: tuple region OOB (need=%d, have=%d)",
                (Int)need_tup, (Int)avail);
    }
    if (data->exp == 0) {
        ErrorQuit("CIB/inc: exp=0 (invalid state)", 0, 0);
    }
    if (data->rank == 0) {
        ErrorQuit("CIB/inc: rank=0 (invalid state)", 0, 0);
    }
#endif

#if defined(CIB_DIAG) && CIB_DIAG>1
  fprintf(stderr,
          "[cib:inc] ctx=%p size=%zu | dim=%zu dim_aligned=%zu exp=%u rank=%zu size(rows)=%zu len=%zu raw_off=%zu tup_off=%zu\n",
          (void*)ADDR_OBJ(ctx), (size_t)INLINE_BYTES(ctx),
          data->dim, data->dim_aligned, (unsigned)data->exp,
          data->rank, data->size, data->len, data->raw_offset, data->tup_offset);
#endif

    uint16_t *c = TUP_PTR(data);
    if ( c[0] ) {
        return false;
    }

    size_t         pos = data->rank;
    const uint16_t max = data->exp - 1;

    while (pos > 0 && c[pos] == max)
    {
        // we do not need to set it to zero, as simple adding one and mod exp will do the job
        add_pos( ctx, pos );

        // po modyfikacji wewnątrz torby odśwież wskaźniki
        // this may be not needed, as we do not fire gc
        // data = DATA_PTR_INLINE(ctx);
        // c    = TUP_PTR(data);

        // reset współrzędnej
        c[pos--] = 0;
    }

    if (pos == 0) {
        c[0] = 1; // overflow
        return false;
    }

    c[pos]++;
    add_pos( ctx, pos );

    return true;
}

Obj FuncCIBVectorSystemContextRawInc(Obj self, Obj ctx)
{
    CHECK_ARGUMENT_CTX(ctx);

#if defined(CIB_DIAG) && CIB_DIAG>1
    fprintf(stderr, "[cib:RawInc] ctxTNUM=%u sz(ctx)=%zu\n",
          (unsigned)TNUM_OBJ(ctx), SIZE_OBJ(ctx));
#endif

    bool result = inc( ctx );
    CHANGED_BAG(ctx);

    return result ? True : False;
}

/* ===== access/display matrices ===== */

Obj FuncCIBVectorSystemContextRawMat(Obj self, Obj ctx, Obj pos)
{
    CHECK_ARGUMENT_CTX(ctx);
    RequireNonnegativeSmallInt( SELF_NAME, pos );

    const size_t npos = INT_INTOBJ(pos);
    cib_data_t  *data = DATA_PTR_INLINE(ctx);

#if defined(CIB_DIAG) && CIB_DIAG
    const size_t avail = INLINE_BYTES(ctx);

    const size_t need_tup = data->tup_offset + (data->rank + 1) * sizeof(uint16_t);
    if (need_tup > avail) {
        ErrorQuit("CIB/inc: tuple region OOB (need=%d, have=%d)",
                (Int)need_tup, (Int)avail);
    }
    if (data->exp == 0) {
        ErrorQuit("CIB/inc: exp=0 (invalid state)", 0, 0);
    }
    if (data->rank == 0) {
        ErrorQuit("CIB/inc: rank=0 (invalid state)", 0, 0);
    }
#endif


    if (npos && npos > data->rank) {
        ErrorQuit("CIBVectorSystemContextRawMat: <pos>=%d out of range", (Int)npos, 0);
    }

    uint16_t    *arr  = NULL;
    const size_t size = data->size;
    const size_t dim  = data->dim;
    const size_t pad  = data->dim_aligned - dim;

    Obj res = NEW_PLIST( T_PLIST, size );
    SET_LEN_PLIST( res, size );
    for (size_t i=1, k=0; i<=size; ++i) {
        Obj row = NEW_PLIST( T_PLIST, dim );
        SET_LEN_PLIST( row, dim );
        for (size_t j=1; j<=dim; ++j) {
            arr = IND_PTR( DATA_PTR_INLINE(ctx), npos );
            SET_ELM_PLIST( row, j, ObjInt_UInt(arr[k++]) );
            CHANGED_BAG(row);
        }
        SET_ELM_PLIST( res, i, row );
        CHANGED_BAG(res);
        k += pad;
    }
    return res;
}

Obj FuncCIBVectorSystemContextRawMatDisplay(Obj self, Obj ctx, Obj pos)
{
    CHECK_ARGUMENT_CTX(ctx);
    RequireNonnegativeSmallInt( SELF_NAME, pos );

    const size_t npos = INT_INTOBJ(pos);
    cib_data_t  *data = DATA_PTR_INLINE(ctx);

    if (npos && npos > data->rank) {
        ErrorQuit("CIBVectorSystemContextRawMat: <pos>=%d out of range", (Int)npos, 0);
    }

    uint16_t    *arr  = IND_PTR(data, npos);
    const size_t size = data->size;
    const size_t dim  = data->dim;
    const size_t pad  = data->dim_aligned - dim;

    for (size_t i=1, k=0; i<=size; ++i) {
        for (size_t j=1; j<=dim; ++j) {
            printf( "%6d", arr[k++] );
        }
        printf("\n");
        k += pad;
    }
    return (Obj)0;
}

/* ===== check if we have a group ===== */

static inline bool equal_rows(const uint16_t* __restrict a, const uint16_t* __restrict b, size_t w)
{
    for (size_t i = 0; i < w; ++i)
        if (a[i] != b[i]) return false;
    return true;
}

static inline bool contains_row(const uint16_t* __restrict base, size_t rows, size_t stride, const uint16_t* __restrict key)
{
    for (size_t r = 0; r < rows; ++r) {
        const uint16_t* row = base + r * stride;
        if (equal_rows(row, key, stride)) return true;
    }
    return false;
}

static inline void neg_row_mod(uint16_t* __restrict dst, const uint16_t* __restrict src, size_t w, uint32_t exp)
{
    for (size_t i = 0; i < w; ++i) {
        uint32_t v = (uint32_t)src[i];
        dst[i] = (uint16_t)(v == 0 ? 0u : (exp - v));
    }
}

static inline void add_rows_mod(uint16_t* __restrict dst, const uint16_t* __restrict a, const uint16_t* __restrict b, size_t w, uint32_t exp)
{
    for (size_t i = 0; i < w; ++i) {
        uint32_t s = (uint32_t)a[i] + (uint32_t)b[i];
        dst[i] = (uint16_t)(s >= exp ? (s - exp) : s);
    }
}

static inline void dbl_row_mod(uint16_t* __restrict dst, const uint16_t* __restrict a, size_t w, uint32_t exp)
{
    for (size_t i = 0; i < w; ++i) {
        uint32_t s = (uint32_t)a[i] + (uint32_t)a[i];
        dst[i] = (uint16_t)(s >= exp ? (s - exp) : s);
    }
}

static inline bool is_group(cib_data_t *data)
{
    const size_t   rows   = data->size;
    const size_t   stride = data->dim_aligned;
    const uint32_t exp    = (uint32_t)data->exp;

    if (rows == 0) return false;

    const uint16_t * __restrict base = RAW_PTR(data);

    // check identity (zero vector)
    bool has_zero = false;
    for (size_t r = 0; r < rows && !has_zero; ++r) {
            const uint16_t *row      = base + r * stride;
            bool            all_zero = true;
            for (size_t j = 0; j < stride; ++j) {
                if (row[j] != 0) { all_zero = false; break; }
            }
            if (all_zero) has_zero = true;
    }
    if (!has_zero) return false;

    // padded temporary buffer for -v, 2v, v+w
    uint16_t tmp[stride];

    for (size_t i = 0; i < rows; ++i) {
        const uint16_t *vi = base + i * stride;

        // opposite
        neg_row_mod(tmp, vi, stride, exp);
        if (!contains_row(base, rows, stride, tmp)) return false;

        // 2*v
        dbl_row_mod(tmp, vi, stride, exp);
        if (!contains_row(base, rows, stride, tmp)) return false;

        for (size_t j = i + 1; j < rows; ++j) {
            const uint16_t *vj = base + j * stride;
            // uniqueness
            if (equal_rows(vi, vj, stride)) return false;
            // closure under addition
            add_rows_mod(tmp, vi, vj, stride, exp);
            if (!contains_row(base, rows, stride, tmp)) return false;
        }
    }
    return true;
}

Obj FuncIsGroupByCIBVectorSystemContextRaw(Obj self, Obj ctx)
{
    CHECK_ARGUMENT_CTX(ctx);
    return is_group( DATA_PTR_INLINE(ctx) ) ? True : False;
}

static StructGVarFunc GVarFuncs [] = {
    GVAR_FUNC(CIBVectorSystemContextCreate,          3, "exp, vs, basis" ),
    GVAR_FUNC(CIBVectorSystemContextRawCoordinates,  1, "ctx" ),
    GVAR_FUNC(CIBVectorSystemContextRawInc,          1, "ctx" ),
    GVAR_FUNC(CIBVectorSystemContextRawSet,          2, "ctx, coordinates"),
    GVAR_FUNC(CIBVectorSystemContextRawMat,          2, "ctx, pos" ),
    GVAR_FUNC(CIBVectorSystemContextRawMatDisplay,   2, "ctx, pos" ),
    GVAR_FUNC(IsGroupByCIBVectorSystemContextRaw,    1, "ctx" ),
    { 0 } /* Finish with an empty entry */
};

static Int InitKernel( StructInitInfo *module )
{
    InitGlobalBag(&CIBVectorSystemContextDataType, "CIBVectorSystemContextDataType");
    ImportGVarFromLibrary("CIBVectorSystemContextDataType", &CIBVectorSystemContextDataType);

    /* set proper VEC_BITS... value */
    VEC_BITS  = detect_vec_bits_runtime();
    VEC_BYTES = VEC_BITS / 8u;

    HEAD_BYTES = find_multiplicity_roof(sizeof(Obj), VEC_BYTES);
    HEAD_SLOTS = HEAD_BYTES / sizeof(Obj);

    /* init filters and functions */
    InitHdlrFuncsFromTable( GVarFuncs );

    /* return success */
    return 0;
}

static Int InitLibrary( StructInitInfo *module )
{
    InitGVarFuncsFromTable( GVarFuncs );
    return 0;
}

static StructInitInfo module = {
    .type       = MODULE_DYNAMIC,
    .name       = "cib",
    .initKernel = InitKernel,
    .initLibrary= InitLibrary,
};

StructInitInfo *Init__Dynamic( void )
{
    return &module;
}