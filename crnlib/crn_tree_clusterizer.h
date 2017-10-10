// File: crn_tree_clusterizer.h
// See Copyright Notice and license at the end of inc/crnlib.h
#pragma once
#include "crn_matrix.h"

namespace crnlib {
template <typename VectorType>
class tree_clusterizer {
 public:
  tree_clusterizer() {}

  struct VectorInfo {
    uint index;
    uint weight;
    float weightedDotProduct;
  };

  void clear() {
    m_hist.clear();
    m_vectors.clear();
    m_codebook.clear();
    m_nodes.clear();
    m_node_index_map.clear();
  }

  void add_training_vec(const VectorType& v, uint weight) {
    m_hist.push_back(std::make_pair(v, weight));
  }

  bool generate_codebook(uint max_size, bool generate_node_index_map = false) {
    if (m_hist.empty())
      return false;

    double ttsum = 0.0f;

    vq_node root;
    root.m_vectors.reserve(static_cast<uint>(m_hist.size()));
    m_vectors.reserve(m_hist.size());

    std::sort(m_hist.begin(), m_hist.end());
    for (uint i = 0; i < m_hist.size(); i++) {
      if (!root.m_vectors.size() || m_hist[i].first != m_vectors.back()) {
        VectorInfo vectorInfo;
        vectorInfo.index = m_vectors.size();
        vectorInfo.weight = m_hist[i].second;
        root.m_vectors.push_back(vectorInfo);
        m_vectors.push_back(m_hist[i].first);
      } else if (root.m_vectors.back().weight > UINT_MAX - m_hist[i].second) {
        root.m_vectors.back().weight = UINT_MAX;
      } else {
        root.m_vectors.back().weight += m_hist[i].second;
      }
    }

    m_weightedVectors.resize(m_vectors.size());
    m_left_children_indices.resize(m_vectors.size());
    m_right_children_indices.resize(m_vectors.size());

    for (uint i = 0; i < m_vectors.size(); i++) {
      const VectorType& v = m_vectors[i];
      const uint weight = root.m_vectors[i].weight;
      m_weightedVectors[i] = v * (float)weight;
      root.m_centroid += m_weightedVectors[i];
      root.m_total_weight += weight;
      root.m_vectors[i].weightedDotProduct = v.dot(v) * weight;
      ttsum += root.m_vectors[i].weightedDotProduct;
    }

    root.m_variance = (float)(ttsum - (root.m_centroid.dot(root.m_centroid) / root.m_total_weight));

    root.m_centroid *= (1.0f / root.m_total_weight);

    m_nodes.clear();
    m_nodes.reserve(max_size * 2 + 1);

    m_nodes.push_back(root);

    // Warning: if this code is NOT compiled with -fno-strict-aliasing, m_nodes.get_ptr() can be NULL here. (Argh!)

    uint total_leaves = 1;

    while (total_leaves < max_size) {
      int worst_node_index = -1;
      float worst_variance = -1.0f;

      for (uint i = 0; i < m_nodes.size(); i++) {
        vq_node& node = m_nodes[i];

        // Skip internal and unsplittable nodes.
        if ((node.m_left != -1) || (node.m_unsplittable))
          continue;

        if (node.m_variance > worst_variance) {
          worst_variance = node.m_variance;
          worst_node_index = i;
        }
      }

      if (worst_variance <= 0.0f)
        break;

      split_node(worst_node_index);
      total_leaves++;
    }

    m_codebook.clear();

    for (uint i = 0; i < m_nodes.size(); i++) {
      vq_node& node = m_nodes[i];
      if (node.m_left != -1) {
        CRNLIB_ASSERT(node.m_right != -1);
        continue;
      }

      CRNLIB_ASSERT((node.m_left == -1) && (node.m_right == -1));

      node.m_codebook_index = m_codebook.size();
      m_codebook.push_back(node.m_centroid);

      if (generate_node_index_map) {
        for (uint j = 0; j < node.m_vectors.size(); j++)
          m_node_index_map.insert(std::make_pair(m_vectors[node.m_vectors[j].index], node.m_codebook_index));
      }
    }

    return true;
  }

  inline uint get_node_index(const VectorType& v) {
    return m_node_index_map.find(v)->second;
  }

  inline uint get_codebook_size() const {
    return m_codebook.size();
  }

  inline const VectorType& get_codebook_entry(uint index) const {
    return m_codebook[index];
  }

  typedef crnlib::vector<VectorType> vector_vec_type;
  inline const vector_vec_type& get_codebook() const {
    return m_codebook;
  }

 private:

  crnlib::vector<std::pair<VectorType, uint> > m_hist;
  crnlib::vector<VectorType> m_vectors;
  crnlib::vector<VectorType> m_weightedVectors;
  crnlib::vector<uint> m_left_children_indices;
  crnlib::vector<uint> m_right_children_indices;
  crnlib::hash_map<VectorType, uint> m_node_index_map;

  struct vq_node {
    vq_node()
        : m_centroid(cClear), m_total_weight(0), m_left(-1), m_right(-1), m_codebook_index(-1), m_unsplittable(false) {}

    VectorType m_centroid;
    uint64 m_total_weight;

    float m_variance;

    crnlib::vector<VectorInfo> m_vectors;

    int m_left;
    int m_right;

    int m_codebook_index;

    bool m_unsplittable;
  };

  typedef crnlib::vector<vq_node> node_vec_type;

  node_vec_type m_nodes;

  vector_vec_type m_codebook;

  void split_node(uint index) {
    vq_node& parent_node = m_nodes[index];

    if (parent_node.m_vectors.size() == 1)
      return;

    VectorType furthest(0);
    double furthest_dist = -1.0f;

    for (uint i = 0; i < parent_node.m_vectors.size(); i++) {
      const VectorType& v = m_vectors[parent_node.m_vectors[i].index];
      double dist = v.squared_distance(parent_node.m_centroid);
      if (dist > furthest_dist) {
        furthest_dist = dist;
        furthest = v;
      }
    }

    VectorType opposite;
    double opposite_dist = -1.0f;

    for (uint i = 0; i < parent_node.m_vectors.size(); i++) {
      const VectorType& v = m_vectors[parent_node.m_vectors[i].index];
      double dist = v.squared_distance(furthest);
      if (dist > opposite_dist) {
        opposite_dist = dist;
        opposite = v;
      }
    }

    VectorType left_child((furthest + parent_node.m_centroid) * .5f);
    VectorType right_child((opposite + parent_node.m_centroid) * .5f);

    if (parent_node.m_vectors.size() > 2) {
      const uint N = VectorType::num_elements;

      matrix<N, N, float> covar;
      covar.clear();

      for (uint i = 0; i < parent_node.m_vectors.size(); i++) {
        const VectorType v = m_vectors[parent_node.m_vectors[i].index] - parent_node.m_centroid;
        const VectorType w = v * (float)parent_node.m_vectors[i].weight;
        for (uint x = 0; x < N; x++) {
          for (uint y = x; y < N; y++)
            covar[x][y] = covar[x][y] + v[x] * w[y];
        }
      }

      float divider = (float)parent_node.m_total_weight;
      for (uint x = 0; x < N; x++) {
        for (uint y = x; y < N; y++) {
          covar[x][y] /= divider;
          covar[y][x] = covar[x][y];
        }
      }

      VectorType axis(1.0f);
      // Starting with an estimate of the principle axis should work better, but doesn't in practice?
      //left_child - right_child);
      //axis.normalize();

      for (uint iter = 0; iter < 10; iter++) {
        VectorType x;

        double max_sum = 0;

        for (uint i = 0; i < N; i++) {
          double sum = 0;

          for (uint j = 0; j < N; j++)
            sum += axis[j] * covar[i][j];

          x[i] = (float)sum;

          max_sum = i ? math::maximum(max_sum, sum) : sum;
        }

        if (max_sum != 0.0f)
          x *= (float)(1.0f / max_sum);

        axis = x;
      }

      axis.normalize();

      VectorType new_left_child(0.0f);
      VectorType new_right_child(0.0f);

      double left_weight = 0.0f;
      double right_weight = 0.0f;

      for (uint i = 0; i < parent_node.m_vectors.size(); i++) {
        const VectorInfo& vectorInfo = parent_node.m_vectors[i];
        const float weight = (float)vectorInfo.weight;
        double t = (m_vectors[vectorInfo.index] - parent_node.m_centroid) * axis;
        if (t < 0.0f) {
          new_left_child += m_weightedVectors[vectorInfo.index];
          left_weight += weight;
        } else {
          new_right_child += m_weightedVectors[vectorInfo.index];
          right_weight += weight;
        }
      }

      if ((left_weight > 0.0f) && (right_weight > 0.0f)) {
        left_child = new_left_child * (float)(1.0f / left_weight);
        right_child = new_right_child * (float)(1.0f / right_weight);
      }
    }

    uint64 left_weight = 0;
    uint64 right_weight = 0;

    uint left_children_indices_count = 0;
    uint right_children_indices_count = 0;

    float prev_total_variance = 1e+10f;

    float left_variance = 0.0f;
    float right_variance = 0.0f;

    // FIXME: Excessive upper limit
    const uint cMaxLoops = 1024;
    for (uint total_loops = 0; total_loops < cMaxLoops; total_loops++) {
      left_children_indices_count = 0;
      right_children_indices_count = 0;

      VectorType new_left_child(cClear);
      VectorType new_right_child(cClear);

      double left_ttsum = 0.0f;
      double right_ttsum = 0.0f;

      left_weight = 0;
      right_weight = 0;

      for (uint i = 0; i < parent_node.m_vectors.size(); i++) {
        const VectorInfo& vectorInfo = parent_node.m_vectors[i];
        double left_dist2 = left_child.squared_distance(m_vectors[vectorInfo.index]);
        double right_dist2 = right_child.squared_distance(m_vectors[vectorInfo.index]);
        if (left_dist2 < right_dist2) {
          m_left_children_indices[left_children_indices_count++] = i;
          new_left_child += m_weightedVectors[vectorInfo.index];
          left_ttsum += vectorInfo.weightedDotProduct;
          left_weight += vectorInfo.weight;
        } else {
          m_right_children_indices[right_children_indices_count++] = i;
          new_right_child += m_weightedVectors[vectorInfo.index];
          right_ttsum += vectorInfo.weightedDotProduct;
          right_weight += vectorInfo.weight;
        }
      }

      if ((!left_weight) || (!right_weight)) {
        parent_node.m_unsplittable = true;
        return;
      }

      left_variance = (float)(left_ttsum - (new_left_child.dot(new_left_child) / left_weight));
      right_variance = (float)(right_ttsum - (new_right_child.dot(new_right_child) / right_weight));

      new_left_child *= (1.0f / left_weight);
      new_right_child *= (1.0f / right_weight);

      left_child = new_left_child;
      right_child = new_right_child;

      float total_variance = left_variance + right_variance;
      if (total_variance < .00001f)
        break;

      if (((prev_total_variance - total_variance) / total_variance) < .00001f)
        break;

      prev_total_variance = total_variance;
    }

    const uint left_child_index = m_nodes.size();
    const uint right_child_index = m_nodes.size() + 1;

    parent_node.m_left = m_nodes.size();
    parent_node.m_right = m_nodes.size() + 1;

    m_nodes.resize(m_nodes.size() + 2);

    // parent_node is invalid now, because m_nodes has been changed

    vq_node& left_child_node = m_nodes[left_child_index];
    vq_node& right_child_node = m_nodes[right_child_index];

    left_child_node.m_centroid = left_child;
    left_child_node.m_total_weight = left_weight;
    left_child_node.m_vectors.resize(left_children_indices_count);
    for (uint i = 0; i < left_children_indices_count; i++)
      left_child_node.m_vectors[i] = parent_node.m_vectors[m_left_children_indices[i]];
    left_child_node.m_variance = left_variance;

    right_child_node.m_centroid = right_child;
    right_child_node.m_total_weight = right_weight;
    right_child_node.m_vectors.resize(right_children_indices_count);
    for (uint i = 0; i < right_children_indices_count; i++)
      right_child_node.m_vectors[i] = parent_node.m_vectors[m_right_children_indices[i]];
    right_child_node.m_variance = right_variance;
  }
};

}  // namespace crnlib
