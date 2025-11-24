#pragma once
#include "transform.h"
#include "ray.h"

namespace gem
{
    struct sphere3f
    {
        float3 c = {};
        float r = 0;

        static sphere3f unit();

        bool GEM_VECTORCALL contains_point(const float3& point) const;

        bool GEM_VECTORCALL intersects_ray(const ray3f& ray, float3* phit, float* thit);

        float3 GEM_VECTORCALL closest_point(const float3& point) const;
    };

    GEM_INLINE sphere3f sphere3f::unit()
    {
        return { { 0.f, 0.f, 0.f }, 1.f };
    }

    GEM_INLINE bool GEM_VECTORCALL sphere3f::contains_point(const float3& point) const
    {
        return length_squared(point - c) <= r * r;
    }

    GEM_INLINE bool GEM_VECTORCALL sphere3f::intersects_ray(const ray3f& ray, float3* p, float* t)
    {
        // NOTE(gortega): https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection.html
        float3 l = c - ray.p;
        float tca = dot(l, ray.v);
        if (tca < 0.f)
        {
            return false;
        }

        float d2 = dot(l, l) - tca * tca;
        if (d2 > r * r)
        {
            return false;
        }
       
        float thc = sqrtf(r * r - d2);
        float t0 = tca - thc;
        float t1 = tca + thc;
        if (t0 > t1)
        {
            float t2 = t0;
            t0 = t1;
            t1 = t2;
        }

        if (t0 < 0)
        {
            t0 = t1;
            if (t0 < 0)
            {
                return false;
            }
        }

        if (t) t[0] = t0;
        if (p) p[0] = ray.p + ray.v * t0;

        return true;
    }

    GEM_INLINE float3 GEM_VECTORCALL sphere3f::closest_point(const float3& point) const
    {
        float3 cp = point - c;
        float l = cp.length();
        if (l <= r)
        {
            return point;
        }

        return c + cp * (1.f / l) * r;
    }
}
