#pragma once
#include "transform.h"
#include "ray.h"
#include <float.h>

namespace gem
{
    struct range2f
    {
        float2 min = { +FLT_MAX, +FLT_MAX };
        float2 max = { -FLT_MAX, -FLT_MAX };

        float2 center() const;

        float2 extent() const;

        float area() const;

        range2f& GEM_VECTORCALL expand(const float2& point);
        
        range2f& GEM_VECTORCALL expand(const range2f& range);

        bool degenerate() const;

        bool GEM_VECTORCALL contains_point(const float2& point) const;
    };

    GEM_INLINE float2 range2f::center() const
    {
        return (max + min) * 0.5f;
    }

    GEM_INLINE float2 range2f::extent() const
    {
        return (max - min) * 0.5f;
    }

    GEM_INLINE float range2f::area() const
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
        float3 min = { +FLT_MAX, +FLT_MAX, +FLT_MAX };
        float3 max = { -FLT_MAX, -FLT_MAX, -FLT_MAX };

        static range3f unit();

        static range3f unit();

        float3 center() const;

        float3 extent() const;
    
        float volume() const;

        range3f& GEM_VECTORCALL expand(const float3& point);

        range3f& GEM_VECTORCALL expand(const range3f& range);

        range3f& GEM_VECTORCALL transform(const transform3f& transform);

        range3f& GEM_VECTORCALL transform(const transform1f& transform);

        bool degenerate() const;

        bool GEM_VECTORCALL contains_point(const float3& point) const;

        bool GEM_VECTORCALL intersects_ray(const ray3f& ray, float3* p_point, float tolerance = 0.00f);
    };

    range3f GEM_VECTORCALL transform_range(const range3f& range, const transform3f& transform);

    range3f GEM_VECTORCALL transform_range(const range3f& range, const transform1f& transform);

    GEM_INLINE range3f range3f::unit()
    {
        return { .min = float3(-0.5f, -0.5f, -0.5f), .max = float3(+0.5f, +0.5f, +0.5f) };
    }

    GEM_INLINE float3 range3f::center() const
    {
        return (max + min) * 0.5f;
    }

    GEM_INLINE float3 range3f::extent() const
    {
        return (max - min) * 0.5f;
    }

    GEM_INLINE float range3f::volume() const
    {
        return (max.x - min.x) * (max.y - min.y) * (max.z * min.z);

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

    GEM_INLINE range3f& GEM_VECTORCALL range3f::transform(const transform3f& transform)
    {
        float3 points[8] =
        {
            { min.x, min.y, min.z },
            { max.x, min.y, min.z },
            { min.x, max.y, min.z },
            { max.x, max.y, min.z },
            { min.x, min.y, max.z },
            { max.x, min.y, max.z },
            { min.x, max.y, max.z },
            { max.x, max.y, max.z },
        };

        min = { FLT_MAX, FLT_MAX, FLT_MAX };
        max = { FLT_MIN, FLT_MIN, FLT_MIN };
        
        expand(transform.transform_point(points[0]));
        expand(transform.transform_point(points[1]));
        expand(transform.transform_point(points[2]));
        expand(transform.transform_point(points[3]));
        expand(transform.transform_point(points[4]));
        expand(transform.transform_point(points[5]));
        expand(transform.transform_point(points[6]));
        expand(transform.transform_point(points[7]));

        return *this;
    }

    GEM_INLINE range3f& GEM_VECTORCALL range3f::transform(const transform1f& transform)
    {
        float3 points[8] =
        {
            { min.x, min.y, min.z },
            { max.x, min.y, min.z },
            { min.x, max.y, min.z },
            { max.x, max.y, min.z },
            { min.x, min.y, max.z },
            { max.x, min.y, max.z },
            { min.x, max.y, max.z },
            { max.x, max.y, max.z },
        };

        min = { FLT_MAX, FLT_MAX, FLT_MAX };
        max = { FLT_MIN, FLT_MIN, FLT_MIN };

        expand(transform.transform_point(points[0]));
        expand(transform.transform_point(points[1]));
        expand(transform.transform_point(points[2]));
        expand(transform.transform_point(points[3]));
        expand(transform.transform_point(points[4]));
        expand(transform.transform_point(points[5]));
        expand(transform.transform_point(points[6]));
        expand(transform.transform_point(points[7]));

        return *this;
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

    GEM_INLINE bool GEM_VECTORCALL range3f::intersects_ray(const ray3f& ray, float3* p_point, float tolerance /*= 0.00f*/)
    {
        https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
        float3 iv =
        {
            1.0f / ray.v.x,
            1.0f / ray.v.y,
            1.0f / ray.v.z
        };

        int sign[3] =
        {
            iv.x < 0.0f,
            iv.y < 0.0f,
            iv.z < 0.0f
        };

        float3 bounds[2] =
        {
            min,
            max
        };

        float tmin, tmax, tymin, tymax, tzmin, tzmax;
        tmin = (bounds[sign[0]].x - ray.p.x) * iv.x;
        tmax = (bounds[1 - sign[0]].x - ray.p.x) * iv.x;
        tymin = (bounds[sign[1]].y - ray.p.y) * iv.y;
        tymax = (bounds[1 - sign[1]].y - ray.p.y) * iv.y;

        if ((tmin > tymax) || (tymin > tmax))
            return false;
        if (tymin > tmin)
            tmin = tymin;
        if (tymax < tmax)
            tmax = tymax;

        tzmin = (bounds[sign[2]].z - ray.p.z) * iv.z;
        tzmax = (bounds[1 - sign[2]].z - ray.p.z) * iv.z;

        if ((tmin > tzmax) || (tzmin > tmax))
            return false;
        if (tzmin > tmin)
            tmin = tzmin;
        if (tzmax < tmax)
            tmax = tzmax;

        if (p_point)
            *p_point = ray.p + tmin * ray.v;

        return true;
    }

    GEM_INLINE range3f GEM_VECTORCALL transform_range(const range3f& range, const transform3f& transform)
    {
        range3f o;
        o.min = { FLT_MAX, FLT_MAX, FLT_MAX };
        o.max = { FLT_MIN, FLT_MIN, FLT_MIN };
        o.expand(transform.transform_point({ range.min.x, range.min.y, range.min.z }));
        o.expand(transform.transform_point({ range.max.x, range.min.y, range.min.z }));
        o.expand(transform.transform_point({ range.min.x, range.max.y, range.min.z }));
        o.expand(transform.transform_point({ range.max.x, range.max.y, range.min.z }));
        o.expand(transform.transform_point({ range.min.x, range.min.y, range.max.z }));
        o.expand(transform.transform_point({ range.max.x, range.min.y, range.max.z }));
        o.expand(transform.transform_point({ range.min.x, range.max.y, range.max.z }));
        o.expand(transform.transform_point({ range.max.x, range.max.y, range.max.z }));
        return o;
    }

    GEM_INLINE range3f GEM_VECTORCALL transform_range(const range3f& range, const transform1f& transform)
    {
        range3f o;
        o.min = { FLT_MAX, FLT_MAX, FLT_MAX };
        o.max = { FLT_MIN, FLT_MIN, FLT_MIN };
        o.expand(transform.transform_point({ range.min.x, range.min.y, range.min.z }));
        o.expand(transform.transform_point({ range.max.x, range.min.y, range.min.z }));
        o.expand(transform.transform_point({ range.min.x, range.max.y, range.min.z }));
        o.expand(transform.transform_point({ range.max.x, range.max.y, range.min.z }));
        o.expand(transform.transform_point({ range.min.x, range.min.y, range.max.z }));
        o.expand(transform.transform_point({ range.max.x, range.min.y, range.max.z }));
        o.expand(transform.transform_point({ range.min.x, range.max.y, range.max.z }));
        o.expand(transform.transform_point({ range.max.x, range.max.y, range.max.z }));
        return o;
    }

    struct range2d
    {
        double2 min = { +DBL_MAX, +DBL_MAX };
        double2 max = { -DBL_MAX, -DBL_MAX };

        double2 center() const;

        double2 extent() const;

        double area() const;

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

    GEM_INLINE double range2d::area() const
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
        double3 min = { +DBL_MAX, +DBL_MAX, +DBL_MAX };
        double3 max = { -DBL_MAX, -DBL_MAX, -DBL_MAX };

        double3 center() const;

        double3 extent() const;

        double volume() const;

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

    GEM_INLINE double range3d::volume() const
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
        int2 min = { INT_MAX, INT_MAX };
        int2 max = { INT_MIN, INT_MIN };

        int2 center() const;

        int2 extent() const;

        int area() const;

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

    GEM_INLINE int range2i::area() const
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
        int3 min = { INT_MAX, INT_MAX, INT_MAX };
        int3 max = { INT_MIN, INT_MIN, INT_MIN };

        int3 center() const;

        int3 extent() const;

        int volume() const;

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

    GEM_INLINE int range3i::volume() const
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