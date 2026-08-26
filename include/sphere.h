#pragma once
#include "transform.h"
#include "ray.h"

namespace gem
{
    struct sphere
    {
        float3 c = {};
        float r = 0;

        static sphere GEM_VECTORCALL transform(const similarity3f& transform, const sphere& sphere);

        static sphere GEM_VECTORCALL transform(const rigid3f& transform, const sphere& sphere);

        static sphere unit();

        bool GEM_VECTORCALL contains(const float3& point) const;

        bool GEM_VECTORCALL intersects(const ray3f& ray, float3* phit, float* thit);

        float3 GEM_VECTORCALL closest(const float3& point) const;
    };

    GEM_INLINE sphere GEM_VECTORCALL sphere::transform(const similarity3f& transform, const sphere& s)
    {
        sphere o;
        o.c = transform.mulp(s.c);
        o.r = transform.s * s.r;
        return o;
    }

    GEM_INLINE sphere GEM_VECTORCALL sphere::transform(const rigid3f& transform, const sphere& s)
    {
        sphere o;
        o.c = transform.mulp(s.c);
        o.r = s.r;
        return o;
    }

    GEM_INLINE sphere sphere::unit()
    {
        return { { 0.f, 0.f, 0.f }, 1.f };
    }

    GEM_INLINE bool GEM_VECTORCALL sphere::contains(const float3& point) const
    {
        float l2 = length_squared(point - c);
        float r2 = r * r;
        return l2 <= r2;
    }

    GEM_INLINE bool GEM_VECTORCALL sphere::intersects(const ray3f& ray, float3* p, float* t)
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

    GEM_INLINE float3 GEM_VECTORCALL sphere::closest(const float3& point) const
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
