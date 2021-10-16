#pragma once
#include "vector.h"

namespace gem
{
    struct range2f
    {
        float2 min, max;

        float2 center() const;

        float2 extent() const;

        float2 area() const;

        bool degenerate() const;

        bool GEM_VECTORCALL contains_point(const float2& point) const;

        range2f& GEM_VECTORCALL expand(const float2& point);

        range2f& GEM_VECTORCALL expand(const range2f& range);
    };

    GEM_INLINE float2 range2f::center() const
    {
        return (max + min) * 0.5f;
    }

    GEM_INLINE float2 range2f::extent() const
    {
        return (max - min) * 0.5f;
    }

    GEM_INLINE float2 range2f::area() const
    {
        return (max.x - min.x) * (max.y - min.y);
    }

    GEM_INLINE bool range2f::degenerate() const
    {
        return min.x > max.x
            || min.y > max.y;
    }

    GEM_INLINE bool GEM_VECTORCALL range2f::contains_point(const float2& point) const
    {
        return min.x <= point.x && point.x <= max.x
            && min.y <= point.y && point.y <= max.y;
    }

    GEM_INLINE range2f& GEM_VECTORCALL range2f::expand(const float2& point)
    {
        if (point.x < min.x) min.x = point.x;
        if (point.y < min.y) min.y = point.y;
        if (point.x > max.x) max.x = point.x;
        if (point.y > max.y) max.y = point.y;
        return *this;
    }

    GEM_INLINE range2f& GEM_VECTORCALL range2f::expand(const range2f& range)
    {
        expand(range.min);
        expand(range.max);
        return *this;
    }

    struct range3f
    {
        typedef float3 float3;

        float3 min, max;

        float3 center() const;

        float3 extent() const;
    
        float3 volume() const;

        bool degenerate() const;

        bool GEM_VECTORCALL contains_point(const float3& point) const;

        range3f& GEM_VECTORCALL expand(const float3& point);

        range3f& GEM_VECTORCALL expand(const range3f& range);
    };

    GEM_INLINE float3 range3f::center() const
    {
        return (max + min) * 0.5f;
    }

    GEM_INLINE float3 range3f::extent() const
    {
        return (max - min) * 0.5f;
    }

    GEM_INLINE float3 range3f::volume() const
    {
        return (max.x - min.x) * (max.y - min.y) * (max.z * min.z);

    }

    GEM_INLINE bool range3f::degenerate() const
    {
        return min.x > max.x 
            || min.y > max.y
            || min.z > min.z;
    }

    GEM_INLINE bool GEM_VECTORCALL range3f::contains_point(const float3& point) const
    {
        return min.x <= point.x && point.x <= max.x
            && min.y <= point.y && point.y <= max.y
            && min.z <= point.z && point.z <= max.z;
    }

    GEM_INLINE range3f& GEM_VECTORCALL range3f::expand(const float3& point)
    {
        if (point.x < min.x) min.x = point.x;
        if (point.y < min.y) min.y = point.y;
        if (point.z < min.z) min.z = point.z;
        if (point.x > max.x) max.x = point.x;
        if (point.y > max.y) max.y = point.y;
        if (point.z > max.z) max.z = point.z;
        return *this;
    }

    GEM_INLINE range3f& GEM_VECTORCALL range3f::expand(const range3f& range)
    {
        expand(range.min);
        expand(range.max);
        return *this;
    }

    struct range2d
    {
        double2 min, max;

        double2 center() const;

        double2 extent() const;

        double2 area() const;

        bool degenerate() const;

        bool GEM_VECTORCALL contains_point(const double2& point) const;

        range2d& GEM_VECTORCALL expand(const double2& point);

        range2d& GEM_VECTORCALL expand(const range2d& range);
    };

    GEM_INLINE double2 range2d::center() const
    {
        return (max + min) * 0.5f;
    }

    GEM_INLINE double2 range2d::extent() const
    {
        return (max - min) * 0.5f;
    }

    GEM_INLINE double2 range2d::area() const
    {
        return (max.x - min.x) * (max.y - min.y);
    }

    GEM_INLINE bool range2d::degenerate() const
    {
        return min.x > max.x
            || min.y > max.y;
    }

    GEM_INLINE bool GEM_VECTORCALL range2d::contains_point(const double2& point) const
    {
        return min.x <= point.x && point.x <= max.x
            && min.y <= point.y && point.y <= max.y;
    }

    GEM_INLINE range2d& GEM_VECTORCALL range2d::expand(const double2& point)
    {
        if (point.x < min.x) min.x = point.x;
        if (point.y < min.y) min.y = point.y;
        if (point.x > max.x) max.x = point.x;
        if (point.y > max.y) max.y = point.y;
        return *this;
    }

    GEM_INLINE range2d& GEM_VECTORCALL range2d::expand(const range2d& range)
    {
        expand(range.min);
        expand(range.max);
        return *this;
    }

    struct range3d
    {
        typedef double3 double3;

        double3 min, max;

        double3 center() const;

        double3 extent() const;

        double3 volume() const;

        bool degenerate() const;

        bool GEM_VECTORCALL contains_point(const double3& point) const;

        range3d& GEM_VECTORCALL expand(const double3& point);

        range3d& GEM_VECTORCALL expand(const range3d& range);
    };

    GEM_INLINE double3 range3d::center() const
    {
        return (max + min) * 0.5f;
    }

    GEM_INLINE double3 range3d::extent() const
    {
        return (max - min) * 0.5f;
    }

    GEM_INLINE double3 range3d::volume() const
    {
        return (max.x - min.x) * (max.y - min.y) * (max.z * min.z);
    }

    GEM_INLINE bool range3d::degenerate() const
    {
        return min.x > max.x
            || min.y > max.y
            || min.z > min.z;
    }

    GEM_INLINE bool GEM_VECTORCALL range3d::contains_point(const double3& point) const
    {
        return min.x <= point.x && point.x <= max.x
            && min.y <= point.y && point.y <= max.y
            && min.z <= point.z && point.z <= max.z;
    }

    GEM_INLINE range3d& GEM_VECTORCALL range3d::expand(const double3& point)
    {
        if (point.x < min.x) min.x = point.x;
        if (point.y < min.y) min.y = point.y;
        if (point.z < min.z) min.z = point.z;
        if (point.x > max.x) max.x = point.x;
        if (point.y > max.y) max.y = point.y;
        if (point.z > max.z) max.z = point.z;
        return *this;
    }

    GEM_INLINE range3d& GEM_VECTORCALL range3d::expand(const range3d& range)
    {
        expand(range.min);
        expand(range.max);
        return *this;
    }

    struct range2i
    {
        int2 min, max;

        int2 center() const;

        int2 extent() const;

        int2 area() const;

        bool degenerate() const;

        bool GEM_VECTORCALL contains_point(const int2& point) const;

        range2i& GEM_VECTORCALL expand(const int2& point);

        range2i& GEM_VECTORCALL expand(const range2i& range);
    };

    GEM_INLINE int2 range2i::center() const
    {
        return (max + min) / 2;
    }

    GEM_INLINE int2 range2i::extent() const
    {
        return (max - min) / 2;
    }

    GEM_INLINE int2 range2i::area() const
    {
        return (max.x - min.x) * (max.y - min.y);
    }

    GEM_INLINE bool range2i::degenerate() const
    {
        return min.x > max.x
            || min.y > max.y;
    }

    GEM_INLINE bool GEM_VECTORCALL range2i::contains_point(const int2& point) const
    {
        return min.x <= point.x && point.x <= max.x
            && min.y <= point.y && point.y <= max.y;
    }

    GEM_INLINE range2i& GEM_VECTORCALL range2i::expand(const int2& point)
    {
        if (point.x < min.x) min.x = point.x;
        if (point.y < min.y) min.y = point.y;
        if (point.x > max.x) max.x = point.x;
        if (point.y > max.y) max.y = point.y;
        return *this;
    }

    GEM_INLINE range2i& GEM_VECTORCALL  range2i::expand(const range2i& range)
    {
        expand(range.min);
        expand(range.max);
        return *this;
    }

    struct range3i
    {
        typedef int3 int3;

        int3 min, max;

        int3 center() const;

        int3 extent() const;

        int3 volume() const;

        bool degenerate() const;

        bool GEM_VECTORCALL contains_point(const int3& point) const;

        range3i& GEM_VECTORCALL expand(const range3i& range);

        range3i& GEM_VECTORCALL expand(const int3& point);
    };

    GEM_INLINE int3 range3i::center() const
    {
        return (max + min) / 2;
    }

    GEM_INLINE int3 range3i::extent() const
    {
        return (max - min) / 2;
    }

    GEM_INLINE int3 range3i::volume() const
    {
        return (max.x - min.x) * (max.y - min.y) * (max.z * min.z);
    }

    GEM_INLINE bool range3i::degenerate() const
    {
        return min.x > max.x
            || min.y > max.y
            || min.z > min.z;
    }

    GEM_INLINE bool GEM_VECTORCALL range3i::contains_point(const int3& point) const
    {
        return min.x <= point.x && point.x <= max.x
            && min.y <= point.y && point.y <= max.y
            && min.z <= point.z && point.z <= max.z;
    }

    GEM_INLINE range3i& GEM_VECTORCALL range3i::expand(const int3& point)
    {
        if (point.x < min.x) min.x = point.x;
        if (point.y < min.y) min.y = point.y;
        if (point.z < min.z) min.z = point.z;
        if (point.x > max.x) max.x = point.x;
        if (point.y > max.y) max.y = point.y;
        if (point.z > max.z) max.z = point.z;
        return *this;
    }

    GEM_INLINE range3i& GEM_VECTORCALL range3i::expand(const range3i& range)
    {
        expand(range.min);
        expand(range.max);
        return *this;
    }
}