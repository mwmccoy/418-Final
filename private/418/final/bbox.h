#pragma once

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <ostream>
#include <vector>

#include "vec3.h"

struct BBox {

    /// Default min is max float value, default max is negative max float value
    BBox() : min(FLT_MAX), max(-FLT_MAX) {}
    /// Set minimum and maximum extent
    explicit BBox(Vec3 min, Vec3 max) : min(min), max(max) {}

    BBox(const BBox &) = default;
    BBox &operator=(const BBox &) = default;
    ~BBox() = default;

    /// Rest min to max float, max to negative max float
    void reset() {
        min = Vec3(FLT_MAX);
        max = Vec3(-FLT_MAX);
    }

    /// Expand bounding box to include point
    void enclose(Vec3 point) {
        min = hmin(min, point);
        max = hmax(max, point);
    }
    void enclose(BBox box) {
        min = hmin(min, box.min);
        max = hmax(max, box.max);
    }

    /// Get center point of box
    Vec3 center() const { return (min + max) * 0.5f; }

    // Check whether box has no volume
    bool empty() const { return min.x > max.x || min.y > max.y || min.z > max.z; }

    /// Get surface area of the box
    float surface_area() const {
        if (empty())
            return 0.0f;
        Vec3 extent = max - min;
        return 2.0f * (extent.x * extent.z + extent.x * extent.y + extent.y * extent.z);
    }

    /// Get the eight corner points of the bounding box
    std::vector<Vec3> corners() const {
        std::vector<Vec3> ret(8);
        ret[0] = Vec3(min.x, min.y, min.z);
        ret[1] = Vec3(max.x, min.y, min.z);
        ret[2] = Vec3(min.x, max.y, min.z);
        ret[3] = Vec3(min.x, min.y, max.z);
        ret[4] = Vec3(max.x, max.y, min.z);
        ret[5] = Vec3(min.x, max.y, max.z);
        ret[6] = Vec3(max.x, min.y, max.z);
        ret[7] = Vec3(max.x, max.y, max.z);
        return ret;
    }

    Vec3 min, max;
};

inline std::ostream &operator<<(std::ostream &out, BBox b) {
    out << "BBox{" << b.min << "," << b.max << "}";
    return out;
}
