// File: crn_zeng.h
// See Copyright Notice and license at the end of inc/crnlib.h

#define CRN_ZENG_USE_SPARSE_ARRAY 1

#if CRN_ZENG_USE_SPARSE_ARRAY
#include "crn_sparse_array.h"
#endif

namespace crnlib {

#if CRN_ZENG_USE_SPARSE_ARRAY
typedef sparse_array<uint, 4> hist_type;
#else
typedef crnlib::vector<uint> hist_type;
#endif

typedef float (*zeng_similarity_func)(uint index_a, uint index_b, void* pContext);
void update_hist(hist_type& hist, int i, int j, int n);
void create_zeng_reorder_table(uint n, hist_type& hist, crnlib::vector<uint>& remap_table, zeng_similarity_func pFunc, void* pContext, float similarity_func_weight);

}  // namespace crnlib
