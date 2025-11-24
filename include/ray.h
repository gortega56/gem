#pragma once
#include "transform.h"

namespace gem
{
    struct ray3f
    {
        float3 p;
        float3 v;
    
        static ray3f GEM_VECTORCALL transform(const transform3f& transform, const ray3f& ray);

        static ray3f GEM_VECTORCALL transform(const transform1f& transform, const ray3f& ray);

        void GEM_VECTORCALL transform(const transform3f& transform);

        void GEM_VECTORCALL transform(const transform1f& transform);

        bool GEM_VECTORCALL intersects_ray(const ray3f& ray, float3* p_point, const float tolerance = 0.001f) const;
    };
    
    GEM_INLINE void GEM_VECTORCALL ray3f::transform(const transform3f& transform)
    {
        p = transform.transform_point(p);
        v = transform.transform_vector(v);
    }

    GEM_INLINE void GEM_VECTORCALL ray3f::transform(const transform1f& transform)
    {
        p = transform.transform_point(p);
        v = transform.transform_vector(v);
    }

    GEM_INLINE bool ray3f::intersects_ray(const ray3f& ray, float3* p_point, const float tolerance /*= 0.001f*/) const
    {
        float3 v0 = v;
        float3 v1 = ray.v;
        float3 p0 = p;
        float3 p1 = ray.p;
        float3 dp = p1 - p0;
        float v00 = length_squared(v0);
        float v11 = length_squared(v1);
        float v01 = dot(v0, v1);
        float det = (v01 * v01) - (v00 * v11);
        if (det < tolerance)
            return false;

        det = 1.0f / det;
        float dpv0 = dot(dp, v0);
        float dpv1 = dot(dp, v1);
        float t = ((v01 * dpv1) - (v11 * dpv0)) * det;

        if (p_point)
            *p_point = p0 + v0 * t;

        return true;
    }

    GEM_INLINE ray3f GEM_VECTORCALL ray3f::transform(const transform3f& transform, const ray3f& ray)
    {
        return { transform.transform_point(ray.p), transform.transform_vector(ray.v) };
    }

    GEM_INLINE ray3f GEM_VECTORCALL ray3f::transform(const transform1f& transform, const ray3f& ray)
    {
        return { transform.transform_point(ray.p), transform.transform_vector(ray.v) };
    }
}