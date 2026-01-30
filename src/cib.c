/*
 * CIB: Cofinite Integral Braces in GAP
 */

#include <gap_all.h>    // GAP headers
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>

static Obj CIBVectorSystemContextDataType;

#define CHECK_ARGUMENT_CTX(ctx) \
    RequireArgumentCondition( SELF_NAME, ctx, IS_DATOBJ(ctx) && TYPE_OBJ(ctx) == CIBVectorSystemContextDataType, \
        "must be a CIB vector system context data object" )

static unsigned VEC_BITS  = 128; // defaults to SSE2
static unsigned VEC_BYTES = 16;  // defaults to SSE2

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

static inline size_t find_multiplicity_roof(const size_t m, const size_t n)
{
    GAP_ASSERT( n>0 );
    return ((m + n - 1) / n) * n;
}

static inline size_t next_pow2(size_t x) {
    if (x <= 2) return 2;
    x--;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
#if SIZE_MAX > 0xFFFFFFFFu
    // 64-bit system: use __builtin_clzl (for unsigned long)
    x |= x >> 32;
#endif
    if (x == SIZE_MAX) {
        // overflow
        return 0;
    }
    return x + 1;
}

static inline size_t compute_aligned_dimension(const size_t dim)
{
    const size_t lanes = VEC_BYTES / sizeof(uint16_t);

    if ( dim > lanes ) {
        return find_multiplicity_roof( dim, lanes );
    }
    return next_pow2(dim);
}

typedef struct {
    size_t    dim;           // dimension of the group
    size_t    dim_aligned;   // aligned dimension
    uint16_t  exp;           // current exponent
    size_t    rank;          // rank of the coboundary group
    size_t    size;          // group order
    size_t    len;           // aligned length of one linearized matrix
    size_t    raw_offset;    // raw data offset
    size_t    tup_offset;    // tuple offset
} cib_data_t;

static inline void dbg_dump_ctx(cib_data_t *d, const char *tag) {
#if defined(GAP_KERNEL_DEBUG) && GAP_KERNEL_DEBUG
    /* Uwaga: fprintf na stderr, żeby nie mieszało się z wyjściem GAP-a */
    fprintf(stderr,
            "[cib:%s] VEC_BYTES=%u  size(rows)=%zu  dim=%zu  dim_aligned=%zu  "
            "len=%zu  rows*stride=%zu  rank=%zu\n",
            tag, VEC_BYTES,
            d->size, d->dim, d->dim_aligned,
            d->len, d->size * d->dim_aligned, d->rank);
#else
    (void)d;
    (void)tag;
#endif
}

static inline uint16_t *RAW_PTR(cib_data_t* data) {
  return (uint16_t*)((uint8_t*)data + data->raw_offset);
}
static inline uint16_t *IND_PTR(cib_data_t *data, const size_t ind) {
	return RAW_PTR(data) + ind * data->len;
}
static inline uint16_t *TUP_PTR(cib_data_t* data) {
    return (uint16_t*)((uint8_t*)data + data->tup_offset);
}

static inline cib_data_t *DATA_PTR(Obj ctx) {
    /* see sha256.c file in GAP src/ directory */
	return (cib_data_t*)(&ADDR_OBJ(ctx)[1]);
}

static inline uint16_t *RAW_PTR_OBJ(Obj ctx) {
    return RAW_PTR( DATA_PTR(ctx) );
}
static inline uint16_t *IND_PTR_OBJ(Obj ctx, const size_t ind) {
    return IND_PTR( DATA_PTR(ctx), ind );
}
static inline uint16_t *TUP_PTR_OBJ(Obj ctx) {
    return TUP_PTR( DATA_PTR(ctx) );
}
static inline uint16_t* ARR_AT(Obj ctx, size_t ind, size_t k) {
  cib_data_t *data = DATA_PTR(ctx);
  GAP_ASSERT(ind <= data->rank);
  GAP_ASSERT(k <= data->len);
  return RAW_PTR(data) + ind * data->len + k;
}

static inline bool list_to_data(Obj list, Obj ctx, const size_t ind)
{
    cib_data_t    *data = DATA_PTR(ctx);
 
    const size_t   size = data->size;
    const size_t    dim = data->dim;
    const uint16_t  exp = data->exp;
    const size_t    pad = data->dim_aligned - dim;

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
            uint16_t *p = ARR_AT( ctx, ind, k );
            *p = (uint16_t)( num>=0 ? num : num + exp);
			k++;
			GAP_ASSERT( k <= data->len );
        }
        k += pad;
		GAP_ASSERT( k <= data->len );
    }
    return true;
}

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
    const size_t dim = LEN_PLIST(row);

    const size_t dim_aligned = compute_aligned_dimension(dim);

    GAP_ASSERT( dim <= dim_aligned );

    // we want to use AVX, hence allocate multiplicity of VEC_BITS bits per one matrix
    const size_t lanes = VEC_BYTES/sizeof(uint16_t);
    const size_t len = find_multiplicity_roof(dim_aligned * size, lanes );
    GAP_ASSERT((len * sizeof(uint16_t)) % VEC_BYTES == 0);
    GAP_ASSERT( len >= size * dim_aligned );
    
    const size_t rank = LEN_PLIST(basis);

    const size_t raw_bytes   = (rank + 1) * len * sizeof(uint16_t); // size of raw data
    const size_t tup_bytes   = (rank + 1) * sizeof(uint16_t);       // size of tuple
    const size_t align_pad   = VEC_BYTES - 1;                       // to align raw data in the GAP Bag
    const size_t total_bytes = sizeof(UInt4)                        // GAP Bag header, after sha256.c in GAP source code
                             + sizeof(cib_data_t)                   // the structure
                             + align_pad + raw_bytes + tup_bytes;   // align plus the data
    
    Obj res           = NewBag(T_DATOBJ, total_bytes);
	SET_TYPE_OBJ(res, CIBVectorSystemContextDataType);

    cib_data_t *data  = DATA_PTR(res);

    data->exp         = INT_INTOBJ(exp);
    data->size        = size;
    data->dim         = dim;
    data->dim_aligned = dim_aligned;
    data->rank        = rank;
    data->len         = len;


    uint8_t *base = (uint8_t*)data + sizeof(cib_data_t); // base address just after the heading (cib_data_t structure)

	uintptr_t raw_addr = ((uintptr_t)base + align_pad) & ~(uintptr_t)(align_pad); // aligned address of the raw data
    data->raw_offset   = (size_t)(raw_addr - (uintptr_t)data);                    // offset of the raw data
    data->tup_offset   = data->raw_offset + raw_bytes;                            // offset of the tuple
	
	dbg_dump_ctx(data, "create1");
    
	memset((void*)((uint8_t*)data + data->raw_offset), 0, raw_bytes);
    memset((void*)((uint8_t*)data + data->tup_offset), 0, tup_bytes);

    if (!list_to_data( vs, res, 0 )) {
        return Fail;
    }

    for (size_t i=1; i<=data->rank; i++) {
        Obj vec = ELM_PLIST(basis, i);
        if (!IS_PLIST(vec)) {
            ErrorQuit("CIBVectorSystemContextCreate: basis elements must be lists", 0, 0);
        }
        if (LEN_PLIST(vec) != data->size) {
            ErrorQuit("CIBVectorSystemContextCreate: basis elements must have the same size as <vs>", 0, 0);
        }
        if (!list_to_data( vec, res, i )) {
			ErrorQuit("CIBVectorSystemContextCreate: invalid basis vector no %d", i, 0);
		}
    }
	dbg_dump_ctx(data, "create2");
    CHANGED_BAG(res);
    return res;
}

Obj FuncCIBVectorSystemContextRawCoordinates(Obj self, Obj raw)
{
    CHECK_ARGUMENT_CTX(raw);

    cib_data_t *data = DATA_PTR(raw);
	uint16_t      *c = TUP_PTR(data);

    if ( c[0] ) {
        return Fail;
    }
    
    Obj vec = NEW_PLIST(T_PLIST, data->rank);
    SET_LEN_PLIST(vec, data->rank);

    for (size_t i=1; i<=data->rank; i++) {
        SET_ELM_PLIST(vec, i, ObjInt_UInt( c[i] ));
        CHANGED_BAG(vec);
    }
    return vec;
}

static inline void add_pos(cib_data_t *data, const size_t ind)
{
    GAP_ASSERT( ind>0 );
    
    uint16_t * __restrict src = IND_PTR(data, ind);
    uint16_t * __restrict dst = RAW_PTR(data);

    const uint32_t exp = data->exp;
    const size_t   len = data->len;
    
    for (size_t i=0; i<len; ++i) {
        uint32_t s = (uint32_t)dst[i] + (uint32_t)src[i];
        
        dst[i] = (uint16_t)(s >= exp ? s - exp : s);
    }
}

static inline bool inc(cib_data_t *data)
{
	uint16_t *c = TUP_PTR(data);
    if ( c[0] ) {
        return false;
    }
    size_t   pos = data->rank;
    uint16_t max = data->exp - 1;

    while (pos > 0 && c[pos] == max) 
    {
        // this in fact cancels all previous additions 
        // of the basis vectors at position pos
        add_pos( data, pos );

        // reset coordinate
        c[pos--] = 0;
    }
    if (pos == 0) {
        c[0] = 1; // mark overflow
        return false;
    }
    c[pos]++;
    add_pos( data, pos );

    return true;
}

Obj FuncCIBVectorSystemContextRawInc(Obj self, Obj raw)
{
    CHECK_ARGUMENT_CTX(raw);
    
    bool result = inc( DATA_PTR(raw) );
    CHANGED_BAG(raw);

    return result ? True : False;
}

Obj FuncCIBVectorSystemContextRawMat(Obj self, Obj raw, Obj pos)
{
    CHECK_ARGUMENT_CTX(raw);

    RequireNonnegativeSmallInt( SELF_NAME, pos );

    cib_data_t *data = DATA_PTR(raw);

    const size_t npos = INT_INTOBJ(pos);
    if (npos && npos >= data->rank) {
        ErrorQuit("CIBVectorSystemContextRawMat: <pos>=%d out of range", npos, 0);
    }

    uint16_t *arr     = IND_PTR(data, npos);

    const size_t size = data->size;
    const size_t dim  = data->dim;
    const size_t pad  = data->dim_aligned - dim;

    Obj res = NEW_PLIST( T_PLIST, size );
    SET_LEN_PLIST( res, size );
    for (size_t i=1, k=0; i<=size; ++i) {
        Obj row = NEW_PLIST( T_PLIST, dim );
        SET_LEN_PLIST( row, dim );
        for (size_t j=1; j<=dim; ++j) {
            SET_ELM_PLIST( row, j, ObjInt_UInt(arr[k++]) );
            CHANGED_BAG(row);
        }
        SET_ELM_PLIST( res, i, row );
        CHANGED_BAG(res);
        k += pad;
    }
    return res;
}

// --- helpers for is_group ---

static inline bool equal_rows(const uint16_t* __restrict a,
                                  const uint16_t* __restrict b,
                                  size_t w)
{
    for (size_t i = 0; i < w; ++i)
        if (a[i] != b[i]) return false;
    return true;
}

static inline bool contains_row(const uint16_t* __restrict base,
                                    size_t rows, size_t stride,
                                    const uint16_t* __restrict key)
{
    for (size_t r = 0; r < rows; ++r) {
        const uint16_t* row = base + r * stride;
        if (equal_rows(row, key, stride)) return true;
    }
    return false;
}

static inline void neg_row_mod(uint16_t* __restrict dst,
                                   const uint16_t* __restrict src,
                                   size_t w, uint32_t exp)
{
    for (size_t i = 0; i < w; ++i) {
        uint32_t v = (uint32_t)src[i];
        dst[i] = (uint16_t)(v == 0 ? 0u : (exp - v));
    }
}

static inline void add_rows_mod(uint16_t* __restrict dst,
                                    const uint16_t* __restrict a,
                                    const uint16_t* __restrict b,
                                    size_t w, uint32_t exp)
{
    for (size_t i = 0; i < w; ++i) {
        uint32_t s = (uint32_t)a[i] + (uint32_t)b[i];
        dst[i] = (uint16_t)(s >= exp ? (s - exp) : s);
    }
}

static inline void dbl_row_mod(uint16_t* __restrict dst,
                                   const uint16_t* __restrict a,
                                   size_t w, uint32_t exp)
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

    dbg_dump_ctx(data, "is_group");
    GAP_ASSERT(rows * stride <= data->len);
    
    // check if we have the identity (zero vector)
    bool has_zero = false;
    for (size_t r = 0; r < rows && !has_zero; ++r) {
        const uint16_t *row = base + r * stride;
        bool all_zero = true;
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

        // check opposite elements
        neg_row_mod(tmp, vi, stride, exp);
        if (!contains_row(base, rows, stride, tmp)) return false;

        // check elements multiplied by 2
        // it is here to check for uniqueness later
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

Obj FuncIsGroupByCIBVectorSystemContextRaw(Obj self, Obj raw)
{
    CHECK_ARGUMENT_CTX(raw);
    return is_group( DATA_PTR(raw) ) ? True : False;
}

// Table of functions to export
static StructGVarFunc GVarFuncs [] = {
    GVAR_FUNC(CIBVectorSystemContextCreate, 3, "exp, vs, basis"),
    GVAR_FUNC(CIBVectorSystemContextRawCoordinates, 1, "ctx"),
    GVAR_FUNC(CIBVectorSystemContextRawInc, 1, "ctx"),
    GVAR_FUNC(CIBVectorSystemContextRawMat, 2, "ctx, pos"),
    GVAR_FUNC(IsGroupByCIBVectorSystemContextRaw, 1, "ctx"),
    { 0 } /* Finish with an empty entry */
};

/****************************************************************************
**
*F  InitKernel( <module> ) . . . . . . . .  initialise kernel data structures
*/
static Int InitKernel( StructInitInfo *module )
{
    ImportGVarFromLibrary("CIBVectorSystemContextDataType", &CIBVectorSystemContextDataType);
    
    /* set proper VEC_BITS value */
    VEC_BITS  = detect_vec_bits_runtime();
    VEC_BYTES = VEC_BITS / 8u;
    
    /* init filters and functions */
    InitHdlrFuncsFromTable( GVarFuncs );

    /* return success */
    return 0;
}

/****************************************************************************
**
*F  InitLibrary( <module> ) . . . . . . .  initialise library data structures
*/
static Int InitLibrary( StructInitInfo *module )
{
    /* init filters and functions */
    InitGVarFuncsFromTable( GVarFuncs );

    /* return success */
    return 0;
}

/****************************************************************************
**
*F  Init__Dynamic() . . . . . . . . . . . . . . . . . table of init functions
*/
static StructInitInfo module = {
    .type = MODULE_DYNAMIC,
    .name = "cib",
    .initKernel = InitKernel,
    .initLibrary = InitLibrary,
};

StructInitInfo *Init__Dynamic( void )
{
    return &module;
}
