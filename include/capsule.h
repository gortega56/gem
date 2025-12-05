#pragma once
#include "transform.h"
#include "ray.h"
#include "sphere.h"
#include "cylinder.h"

namespace gem
{
    struct capsule3f
    {
        float3 min = {};
        float3 max = {};
        float r = 0;

        static capsule3f GEM_VECTORCALL transform(const transform3f& transform, const capsule3f& capsule);

        static capsule3f GEM_VECTORCALL transform(const transform1f& transform, const capsule3f& capsule);

        static capsule3f unit();

        float3 center() const;

        float3 axis() const;

        void GEM_VECTORCALL transform(const transform3f& transform);

        void GEM_VECTORCALL transform(const transform1f& transform);

        bool GEM_VECTORCALL contains_point(const float3& point) const;

        bool GEM_VECTORCALL intersects_ray(const ray3f& ray, float3* phit, float* thit, float tolerance = 0.01f);

        float3 GEM_VECTORCALL closest_point(const float3& point) const;

        float3 GEM_VECTORCALL support(const float3& d) const;
    };

    GEM_INLINE capsule3f GEM_VECTORCALL capsule3f::transform(const transform3f& transform, const capsule3f& capsule)
    {
        // NOTE(gortega): Ignore scale for radius since nonuniform scale would no longer be min capsule
        capsule3f o;
        o.min = transform.transform_point(capsule.min);
        o.max = transform.transform_point(capsule.max);
        o.r = capsule.r;
        return o;
    }

    GEM_INLINE capsule3f GEM_VECTORCALL capsule3f::transform(const transform1f& transform, const capsule3f& capsule)
    {
        capsule3f o;
        o.min = transform.transform_point(capsule.min);
        o.max = transform.transform_point(capsule.max);
        o.r = capsule.r * transform.s;
        return o;
    }

    GEM_INLINE capsule3f capsule3f::unit()
    {
        return {
            { 0, 0, -0.5f },
            { 0, 0, +0.5f },
            1.f
        };
    }

    GEM_INLINE float3 capsule3f::center() const
    {
        return (min + max) * 0.5f;
    }

    GEM_INLINE float3 capsule3f::axis() const
    {
        return max - min;
    }

    GEM_INLINE void GEM_VECTORCALL capsule3f::transform(const transform3f& transform)
    {
        min = transform.transform_point(min);
        max = transform.transform_point(max);
        // NOTE(gortega): Ignore scale for radius since nonuniform scale would no longer be min capsule
    }

    GEM_INLINE void GEM_VECTORCALL capsule3f::transform(const transform1f& transform)
    {
        min = transform.transform_point(min);
        max = transform.transform_point(max);
        r = r * transform.s;
    }

    GEM_INLINE bool GEM_VECTORCALL capsule3f::contains_point(const float3& point) const
    {
        float3 c = point;
        float3 ab = max - min;
        float3 ac = c - min;
        float t = dot(ac, ab) / length_squared(ab);
        if (t < 0.f) t = 0.f;
        if (t > 1.f) t = 1.f;
        float3 d = lerp(min, max, t);
        float l = length_squared(d - c);
        return l <= (r * r);
    }

    GEM_INLINE bool GEM_VECTORCALL capsule3f::intersects_ray(const ray3f& ray, float3* p, float* t, float tolerance /*= 0.01f*/)
    {
        sphere3f c0 = { min, r };
        sphere3f c1 = { max, r };
        cylinder3f c2 = { min, max, r };
        float3 p0, p1, p2;
        float t0, t1, t2;
        bool b0, b1, b2;
        b0 = c0.intersects_ray(ray, &p0, &t0);
        b1 = c1.intersects_ray(ray, &p1, &t1);
        b2 = c2.intersects_ray(ray, &p2, &t2, tolerance);
        float tmin = (t0 < t1) ? ((t0 < t2) ? t0 : t2) : ((t1 < t2) ? t1 : t2);
        if (t) t[0] = tmin;
        if (p) p[0] = ray.p * ray.v * tmin;
        return b0 || b1 || b2;
    }

    GEM_INLINE float3 GEM_VECTORCALL capsule3f::closest_point(const float3& point) const
    {
        float3 c = point;
        float3 ab = max - min;
        float3 ac = c - min;
        float t = dot(ac, ab) / length_squared(ab);
        if (t < 0.f)
        {
            return sphere3f{ min, r }.closest_point(point);
        }
        else if (t > 1)
        {
            return sphere3f{ max, r }.closest_point(point);
        }

        float3 ac_para_ab = lerp(min, max, t);
        float3 ac_perb_ab = ac - ac_para_ab;
        float l = ac_perb_ab.length();
        ac_perb_ab *= (1.f / l);
        if (l > r) l = r;
        return ac_para_ab + ac_perb_ab * l;
    }

    GEM_INLINE float3 GEM_VECTORCALL capsule3f::support(const float3& d) const
    {
        float dab = dot(d, axis());
        float3 o = (dab < 0.f)
            ? min
            : max;
        return o + gem::normalize(d) * r;
    }
}