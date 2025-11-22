#pragma once
#include "transform.h"
#include "ray.h"
#include "sphere.h"
#include "cylinder.h"

namespace gem
{
    struct capsule3f
    {
        float3 a = {};
        float3 b = {};
        float r = 0;

        static capsule3f unit();

        float3 center() const;

        capsule3f& GEM_VECTORCALL transform(const transform3f& transform);

        capsule3f& GEM_VECTORCALL transform(const transform1f& transform);

        bool GEM_VECTORCALL contains_point(const float3& point) const;

        bool GEM_VECTORCALL intersects_ray(const ray3f& ray, float3* phit, float* thit, float tolerance = 0.01f);

        float3 GEM_VECTORCALL closest_point(const float3& point) const;
    };

    GEM_INLINE capsule3f capsule3f::unit()
    {
        return {
            { 0, 0, -1 },
            { 0, 0, +1 },
            1.f
        };
    }

    GEM_INLINE float3 capsule3f::center() const
    {
        return (a + b) * 0.5f;
    }

    GEM_INLINE capsule3f& GEM_VECTORCALL capsule3f::transform(const transform3f& transform)
    {
        a = transform.transform_point(a);
        b = transform.transform_point(b);
        // NOTE(gortega): Ignore scale for radius since nonuniform scale would no longer be a capsule
    }

    GEM_INLINE capsule3f& GEM_VECTORCALL capsule3f::transform(const transform1f& transform)
    {
        a = transform.transform_point(a);
        b = transform.transform_point(b);
        r = r * transform.s;
    }

    GEM_INLINE bool GEM_VECTORCALL capsule3f::contains_point(const float3& point) const
    {
        float3 c = point;
        float3 ab = b - a;
        float3 ac = c - a;
        float t = dot(ac, ab) / length_squared(ab);
        if (t < 0.f) t = 0.f;
        if (t > 1.f) t = 1.f;
        float3 d = lerp(a, b, t);
        float l = length_squared(d - c);
        return l <= (r * r);
    }

    GEM_INLINE bool GEM_VECTORCALL capsule3f::intersects_ray(const ray3f& ray, float3* p, float* t, float tolerance /*= 0.01f*/)
    {
        sphere3f c0 = { a, r };
        sphere3f c1 = { b, r };
        cylinder3f c2 = { a, b, r };
        float3 p0, p1, p2;
        float t0, t1, t2;
        bool b0, b1, b2;
        b0 = c0.intersects_ray(ray, &p0, &t0);
        b1 = c1.intersects_ray(ray, &p1, &t1);
        b2 = c2.intersects_ray(ray, &p2, &t2, tolerance);
        t[0] = (t0 < t1) ? ((t0 < t2) ? t0 : t2) : ((t1 < t2) ? t1 : t2);
        p[0] = ray.p * ray.v * t[0];
        return b0 || b1 || b2;
    }

    GEM_INLINE float3 GEM_VECTORCALL capsule3f::closest_point(const float3& point) const
    {
        float3 c = point;
        float3 ab = b - a;
        float3 ac = c - a;
        float t = dot(ac, ab) / length_squared(ab);
        if (t < 0.f)
        {
            return sphere3f{ a, r }.closest_point(point);
        }
        else if (t > 1)
        {
            return sphere3f{ b, r }.closest_point(point);
        }

        float3 ac_para_ab = lerp(a, b, t);
        float3 ac_perb_ab = ac - ac_para_ab;
        float l = ac_perb_ab.length();
        ac_perb_ab *= (1.f / l);
        if (l > r) l = r;
        return ac_para_ab + ac_perb_ab * l;
    }
}