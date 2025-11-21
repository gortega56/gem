#pragma once
#include "transform.h"
#include "ray.h"
#include <float.h>

namespace gem
{
    struct cylinder3f
    {
        float3 a = {};
        float3 b = {};
        float r = 0;

        float3 center() const;

        cylinder3f& GEM_VECTORCALL transform(const transform3f& transform);

        cylinder3f& GEM_VECTORCALL transform(const transform1f& transform);

        bool GEM_VECTORCALL contains_point(const float3& point) const;

        bool GEM_VECTORCALL intersects_ray(const ray3f& ray, float3* phit, float* thit, float tolerance = 0.01f);

        float3 GEM_VECTORCALL closest_point(const float3& point) const;
    };

    GEM_INLINE float3 cylinder3f::center() const
    {
        return (a + b) * 0.5f;
    }

    GEM_INLINE cylinder3f& GEM_VECTORCALL cylinder3f::transform(const transform3f& transform)
    {
        a = transform.transform_point(a);
        b = transform.transform_point(b);

        // NOTE(gortega): Ignore scale for radius since nonuniform scale would no longer be a capsule
    }

    GEM_INLINE cylinder3f& GEM_VECTORCALL cylinder3f::transform(const transform1f& transform)
    {
        a = transform.transform_point(a);
        b = transform.transform_point(b);
        r = r * transform.s;
    }

    GEM_INLINE bool GEM_VECTORCALL cylinder3f::contains_point(const float3& point) const
    {
        float3 c = point;
        float3 ab = b - a;
        float3 ac = c - a;
        float t = dot(ac, ab) / length_squared(ab);
        bool result = false;
        if (0.f <= t && t <= 1)
        {
            float3 d = lerp(a, b, t);
            bool l = length_squared(d - c);
            result = l <= (r * r);
        }

        return result;
    }

    GEM_INLINE bool GEM_VECTORCALL cylinder3f::intersects_ray(const ray3f& ray, float3* p, float* t, float tolerance /*= 0.01f*/)
    {
        // NOTE(gortega): RTCD p197
        float3 d = b - a;
        float3 m = ray.p - a;
        float3 n = ray.v;
        float md = dot(m, d);
        float nd = dot(n, d);
        float dd = dot(d, d);

        // Test if ray is fully outside each end of the cylinder
        if (md < 0.f && md + nd < 0.f) return false;
        if (md > dd && md + nd > dd) return false;

        float nn = dot(n, n);
        float mn = dot(m, n);
        float a = dd * nn - nd * nd;
        float k = dot(m, m) - r * r;
        float c = dd * k - md * md;

        // Test if ray is parallel
        if (fabsf(a) < tolerance)
        {
            // Test if outide
            if (c > 0.f)
            {
                return false;
            }

            if (md < 0.f)
            {
                t[0] = -mn / nn;
            }
            else if (md > dd)
            {
                t[0] = (nd - mn) / nn;
            }
            else
            {
                t[0] = 0.f;
            }

            return true;
        }

        float b = dd * mn - nd * md;
        float discr = b * b - a * c;
        if (discr < 0.f)
        {
            // No real roots
            return false;
        }

        t[0] = (-b - sqrtf(discr)) / a;
        if (t[0] < 0.f)
        {
            return false;
        }

        if (md + t[0] * nd < 0.f)
        {
            if (nd <= 0.f)
            {
                return false;
            }

            t[0] = -md / nd;

            return k + 2 * t[0] * (mn + t[0] * nn) <= 0.f;
        }
        else if (md + t[0] * nd > dd)
        {
            if (nd >= 0.f)
            {
                return false;
            }

            t[0] = (dd - md) / nd;

            return k + dd - 2 * md + t[0] * (2 * (mn - nd) + t[0] * nn) <= 0.f;
        }

        p[0] = ray.p * ray.v * t[0];

        return true;
    }

    GEM_INLINE float3 GEM_VECTORCALL cylinder3f::closest_point(const float3& point) const
    {
        float3 c = point;
        float3 ab = b - a;
        float3 ac = c - a;
        float t = dot(ac, ab) / length_squared(ab);
        if (t < 0.f) t = 0.f;
        if (t > 1.f) t = 1.f;
        float3 ac_para_ab = lerp(a, b, t);
        float3 ac_perb_ab = ac - ac_para_ab;
        float l = ac_perb_ab.length();
        ac_perb_ab *= (1.f / l);
        if (l > r) l = r;
        return ac_para_ab + ac_perb_ab * l;
    }
}