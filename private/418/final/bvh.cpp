#include "bvh.h"

#include <stack>
#include <cfloat>
#include <iostream>

namespace PT {

/*inline float SAH(auto start, auto mid, auto end, float big_SA) {
    BBox left;
    BBox right;

    for (auto it = start; it != mid; it++) {
        left.enclose(it->bbox());
    }

    for (auto it = mid; it != end; it++) {
        right.enclose(it->bbox());
    }

    float left_cost = left.surface_area()*((float)std::distance(start, mid));
    float right_cost = right.surface_area()*((float)std::distance(mid, end));

    return (left_cost + right_cost)/big_SA;
}*/


template <typename Primitive>
inline void BVH<Primitive>::fill_node(
                std::vector<Primitive>::iterator start,
                std::vector<Primitive>::iterator end,
                size_t curr_node, size_t max_leaf_size) {

    size_t num_buckets = 16;
    // Setting start and size of node
    nodes[curr_node].start = std::distance(primitives.begin(), start);
    nodes[curr_node].size = std::distance(start, end);
    //std::cout << nodes[curr_node].size << std::endl;

    // Create bounding box for all prims in node
    BBox big_box;
    for (auto it = start; it != end; it++) {
        big_box.enclose(it->bbox());
    }
    nodes[curr_node].bbox = big_box;

    // Set node as leaf and return if prim count is less than leaf size
    // std::cout << max_leaf_size << std::endl;
    if ((size_t)std::distance(start, end) <= max_leaf_size) {
        nodes[curr_node].l = 0;
        nodes[curr_node].r = 0;
        return;
    }

    Vec3 box_min = big_box.min;
    Vec3 box_max = big_box.max;


    // Find best split on X
    float best_x = 0.0f;
    float best_xsah = FLT_MAX;

    std::vector<BBox> x_bboxes;
    std::vector<size_t> x_counts;

    for (size_t i = 0; i < num_buckets; i++) {
        x_bboxes.push_back(BBox());
        x_counts.push_back(0);
    }


    for (auto it = start; it != end; it++) {
        size_t idx = std::floor(((float)num_buckets)*(it->bbox().center()[0]
                        - big_box.min[0])/(big_box.max[0]-big_box.min[0]));

        idx = 0;
        for (size_t i = 1; i <= num_buckets; i++) {
            float curr_b = box_min[0] + i*(box_max[0] - box_min[0])/num_buckets;
            if (it->bbox().center()[0] >= curr_b) idx++;
        }

        if (idx == num_buckets) idx--;
        x_bboxes[idx].enclose(it->bbox());
        x_counts[idx]++;
    }

    for (size_t i = 1; i < num_buckets; i++) {
        BBox lbox = BBox();
        BBox rbox = BBox();
        size_t lcount = 0;
        size_t rcount = 0;

        for (size_t j = 0; j != i; j++) {
            lbox.enclose(x_bboxes[j]);
            lcount = lcount + x_counts[j];
        }
        for (size_t k = i; k != num_buckets; k++) {
            rbox.enclose(x_bboxes[k]);
            rcount = rcount + x_counts[k];
        }

        float curr_sah = (rcount*rbox.surface_area() + lcount*lbox.surface_area())/big_box.surface_area();
        if (curr_sah <= best_xsah) {
            best_xsah = curr_sah;
            best_x = box_min[0] + i*(box_max[0] - box_min[0])/num_buckets;
        }
    }
    /*
    for (size_t i = 1; i < num_buckets; i++) {
        float curr_bound = box_min[0] + i*(box_max[0] - box_min[0])/num_buckets;
        auto mid_it = std::partition(start, end, [curr_bound](const Primitive &p)
        {
            return p.bbox().center()[0] < curr_bound;
        });
        float curr_sah = SAH(start, mid_it, end, big_box.surface_area());
        if (curr_sah <= best_xsah) {
            best_xsah = curr_sah;
            best_x = curr_bound;
        }
    }*/

    // Find best split on Y
    float best_y = 0.0f;
    float best_ysah = FLT_MAX;
    size_t mid_dy = 0;
    (void) mid_dy;


    std::vector<BBox> y_bboxes;
    std::vector<size_t> y_counts;

    for (size_t i = 0; i < num_buckets; i++) {
        y_bboxes.push_back(BBox());
        y_counts.push_back(0);
    }

    for (auto it = start; it != end; it++) {
        size_t idx = 0;
        for (size_t i = 1; i <= num_buckets; i++) {
            float curr_b = box_min[1] + i*(box_max[1] - box_min[1])/num_buckets;
            if (it->bbox().center()[1] >= curr_b) idx++;
        }

        if (idx == num_buckets) idx--;
        y_bboxes[idx].enclose(it->bbox());
        y_counts[idx]++;
    }

    for (size_t i = 1; i < num_buckets; i++) {
        BBox lbox = BBox();
        BBox rbox = BBox();
        size_t lcount = 0;
        size_t rcount = 0;

        for (size_t j = 0; j != i; j++) {
            lbox.enclose(y_bboxes[j]);
            lcount = lcount + y_counts[j];
        }
        for (size_t k = i; k != num_buckets; k++) {
            rbox.enclose(y_bboxes[k]);
            rcount = rcount + y_counts[k];
        }

        float curr_sah = (rcount*rbox.surface_area() + lcount*lbox.surface_area())/big_box.surface_area();
        if (curr_sah <= best_ysah) {
            best_ysah = curr_sah;
            best_y = box_min[1] + i*(box_max[1] - box_min[1])/num_buckets;
        }
    }
    /*
    for (size_t i = 1; i < num_buckets; i++) {
        float curr_bound = box_min[1] + i*(box_max[1] - box_min[1])/num_buckets;
        auto mid_it = std::partition(start, end, [curr_bound](const Primitive &p)
        {
            return p.bbox().center()[1] < curr_bound;
        });
        float curr_sah = SAH(start, mid_it, end, big_box.surface_area());
        if (curr_sah <= best_ysah) {
            best_ysah = curr_sah;
            best_y = curr_bound;
            mid_dy = std::distance(start, mid_it);
        }
    }*/

    // Find best split on Z
    float best_z = 0.0f;
    float best_zsah = FLT_MAX;
    size_t mid_dz = 0;
    (void) mid_dz;

    std::vector<BBox> z_bboxes;
    std::vector<size_t> z_counts;

    for (size_t i = 0; i < num_buckets; i++) {
        z_bboxes.push_back(BBox());
        z_counts.push_back(0);
    }

    for (auto it = start; it != end; it++) {
        size_t idx = 0;
        for (size_t i = 1; i <= num_buckets; i++) {
            float curr_b = box_min[2] + i*(box_max[2] - box_min[2])/num_buckets;
            if (it->bbox().center()[2] >= curr_b) idx++;
        }

        if (idx == num_buckets) idx--;
        z_bboxes[idx].enclose(it->bbox());
        z_counts[idx]++;
    }

    for (size_t i = 1; i < num_buckets; i++) {
        BBox lbox = BBox();
        BBox rbox = BBox();
        size_t lcount = 0;
        size_t rcount = 0;

        for (size_t j = 0; j != i; j++) {
            lbox.enclose(z_bboxes[j]);
            lcount = lcount + z_counts[j];
        }
        for (size_t k = i; k != num_buckets; k++) {
            rbox.enclose(z_bboxes[k]);
            rcount = rcount + z_counts[k];
        }

        float curr_sah = (rcount*rbox.surface_area() + lcount*lbox.surface_area())/big_box.surface_area();
        if (curr_sah <= best_zsah) {
            best_zsah = curr_sah;
            best_z = box_min[2] + i*(box_max[2] - box_min[2])/num_buckets;
        }
    }
    /*
    for (size_t i = 1; i < num_buckets; i++) {
        float curr_bound = box_min[2] + i*(box_max[2] - box_min[2])/num_buckets;
        auto mid_it = std::partition(start, end, [curr_bound](const Primitive &p)
        {
            return p.bbox().center()[2] < curr_bound;
        });
        float curr_sah = SAH(start, mid_it, end, big_box.surface_area());
        if (curr_sah <= best_zsah) {
            best_zsah = curr_sah;
            best_z = curr_bound;
            mid_dz = std::distance(start, mid_it);
        }
    }*/

    float best_sah = std::min(best_zsah, std::min(best_ysah, best_xsah));

    if (best_sah == best_xsah) {
        auto mid_it = std::partition(start, end, [best_x](const Primitive &p)
        {
            return p.bbox().center()[0] < best_x;
        });

        size_t d = std::distance(start, mid_it);

        if ((d == 0) || (d == (size_t)std::distance(start, end))) {
            std::sort(start, end, [](const Primitive &p1, const Primitive &p2)
            {
                return p1.bbox().center()[0] < p1.bbox().center()[0];
            });
            mid_it = start + std::distance(start, end)/2;
        }

        size_t left = new_node();
        size_t right = new_node();
        nodes[curr_node].l = left;
        nodes[curr_node].r = right;

        fill_node(start, mid_it, left, max_leaf_size);
        fill_node(mid_it, end, right, max_leaf_size);

    } else if (best_sah == best_ysah) {
        auto mid_it = std::partition(start, end, [best_y](const Primitive &p)
        {
            return p.bbox().center()[1] < best_y;
        });

        size_t d = std::distance(start, mid_it);

        if ((d == 0) || (d == (size_t)std::distance(start, end))) {
            std::sort(start, end, [](const Primitive &p1, const Primitive &p2)
            {
                return p1.bbox().center()[1] < p1.bbox().center()[1];
            });
            mid_it = start + std::distance(start, end)/2;
        }

        size_t left = new_node();
        size_t right = new_node();
        nodes[curr_node].l = left;
        nodes[curr_node].r = right;

        fill_node(start, mid_it, left, max_leaf_size);
        fill_node(mid_it, end, right, max_leaf_size);

    } else {
        auto mid_it = std::partition(start, end, [best_z](const Primitive &p)
        {
            return p.bbox().center()[2] < best_z;
        });

        size_t d = std::distance(start, mid_it);

        if ((d == 0) || (d == (size_t)std::distance(start, end))) {
            std::sort(start, end, [](const Primitive &p1, const Primitive &p2)
            {
                return p1.bbox().center()[2] < p1.bbox().center()[2];
            });
            mid_it = start + std::distance(start, end)/2;
        }

        size_t left = new_node();
        size_t right = new_node();
        nodes[curr_node].l = left;
        nodes[curr_node].r = right;

        fill_node(start, mid_it, left, max_leaf_size);
        fill_node(mid_it, end, right, max_leaf_size);

    }

}


template <typename Primitive>
void BVH<Primitive>::build(std::vector<Primitive> &&prims, size_t max_leaf_size) {

    // NOTE (PathTracer):
    // This BVH is parameterized on the type of the primitive it contains. This allows
    // us to build a BVH over any type that defines a certain interface. Specifically,
    // we use this to both build a BVH over triangles within each Tri_Mesh, and over
    // a variety of Objects (which might be Tri_Meshes, Spheres, etc.) in Pathtracer.
    //
    // The Primitive interface must implement these two functions:
    //      BBox bbox() const;
    //      Trace hit(const Ray& ray) const;
    // Hence, you may call bbox() and hit() on any value of type Primitive.
    //
    // Finally, also note that while a BVH is a tree structure, our BVH nodes don't
    // contain pointers to children, but rather indicies. This is because instead
    // of allocating each node individually, the BVH class contains a vector that
    // holds all of the nodes. Hence, to get the child of a node, you have to
    // look up the child index in this vector (e.g. nodes[node.l]). Similarly,
    // to create a new node, don't allocate one yourself - use BVH::new_node, which
    // returns the index of a newly added node.

    nodes.clear();
    primitives = std::move(prims);

    // TODO (PathTracer): Task 3
    // Construct a BVH from the given vector of primitives and maximum leaf
    // size configuration. The starter code builds a BVH with a
    // single leaf node (which is also the root) that encloses all the
    // primitives.

    size_t head = new_node();
    fill_node(primitives.begin(), primitives.end(), head, max_leaf_size);
}


template <typename Primitive>
BVH<Primitive>::BVH(std::vector<Primitive> &&prims, size_t max_leaf_size) {
    build(std::move(prims), max_leaf_size);
}

template <typename Primitive> bool BVH<Primitive>::Node::is_leaf() const {
    return l == 0 && r == 0;
}

template <typename Primitive>
size_t BVH<Primitive>::new_node(BBox box, size_t start, size_t size, size_t l, size_t r) {
    Node n;
    n.bbox = box;
    n.start = start;
    n.size = size;
    n.l = l;
    n.r = r;
    nodes.push_back(n);
    return nodes.size() - 1;
}

template <typename Primitive> BBox BVH<Primitive>::bbox() const { return nodes[0].bbox; }

template <typename Primitive> std::vector<Primitive> BVH<Primitive>::destructure() {
    nodes.clear();
    return std::move(primitives);
}

template <typename Primitive> void BVH<Primitive>::clear() {
    nodes.clear();
    primitives.clear();
}

} // namespace PT
