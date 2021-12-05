#pragma once

namespace PT {

template <typename Primitive> class BVH {
public:
    BVH() = default;
    BVH(std::vector<Primitive> &&primitives, size_t max_leaf_size = 1);
    void build(std::vector<Primitive> &&primitives, size_t max_leaf_size = 1);
    void fill_node(typename std::vector<Primitive>::iterator start, typename std::vector<Primitive>::iterator end, size_t curr_node, size_t max_leaf_size);


    BBox bbox() const;


    std::vector<Primitive> destructure();
    void clear();

private:
    class Node {
        BBox bbox;
        size_t start, size, l, r;

        bool is_leaf() const;
        friend class BVH<Primitive>;
    };
    size_t new_node(BBox box = {}, size_t start = 0, size_t size = 0, size_t l = 0, size_t r = 0);

    std::vector<Node> nodes;
    std::vector<Primitive> primitives;

    size_t root_idx = 0;
};

} // namespace PT
