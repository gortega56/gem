#pragma once
#include "line.h"
#include "ray.h"

namespace gem
{
    struct plane4f
    {
        float x, y, z, w;

        float distance_from_origin() const;

        float3 normal() const;

        plane4f& normalize();
        
        plane4f& transform(const transform3f& transform);

        plane4f& transform(const transform1f& transform);

        bool GEM_VECTORCALL intersects_plane(const plane4f& plane, line3f* p_line, float tolerance = 0.001f) const;

        bool GEM_VECTORCALL intersects_planes(const plane4f& p0, const plane4f& p1, float3* p_point, float tolerance = 0.001f) const;

        bool GEM_VECTORCALL intersects_ray(const ray3f& ray, float3* p_point, const float tolerance = 0.001f);
    };

    GEM_INLINE float plane4f::distance_from_origin() const
    {
        return fabsf(w) / length(x, y, z);
    }

    GEM_INLINE float3 plane4f::normal() const
    {
        return { x, y, z };
    }

    GEM_INLINE plane4f& plane4f::normalize()
    {
        float il = 1.0f / sqrtf((x * x) + (y * y) + (z * z));
        x *= il;
        y *= il;
        z *= il;
        w *= il;
        return *this;
    }

    GEM_INLINE plane4f& GEM_VECTORCALL plane4f::transform(const transform3f& transform)
    {
        float4x4 m = transform.matrix4x4().inverse().transpose();
        float4 q = { x, y, z, w };
        float4 r = q * m;
        x = r.x;
        y = r.y;
        z = r.z;
        w = r.w;
        return *this;
    }

    GEM_INLINE plane4f& GEM_VECTORCALL plane4f::transform(const transform1f& transform)
    {
        float4x4 m = transform.matrix4x4().inverse().transpose();
        float4 q = { x, y, z, w };
        float4 r = q * m;
        x = r.x;
        y = r.y;
        z = r.z;
        w = r.w;
        return *this;
    }

    GEM_INLINE bool GEM_VECTORCALL plane4f::intersects_plane(const plane4f& plane, line3f* p_line, float tolerance /*= 0.001f*/) const
    {
        float3 n0 = normal();
        float3 n1 = plane.normal();
        float3 v = cross(n0, n1);
        float det = length_squared(v);
        if (det < tolerance)
            return false;

        if (p_line)
            *p_line = { v, (w * n1) - (plane.w * n0) };

        return true;
    }

    GEM_INLINE bool GEM_VECTORCALL plane4f::intersects_planes(const plane4f& p1, const plane4f& p2, float3* p_point, float tolerance /*= 0.001f*/) const
    {
        float3 n0 = normal();
        float3 n1 = p1.normal();
        float3 n2 = p2.normal();
        float3 n0xn1 = cross(n0, n1);
        float det = dot(n0xn1, n2);
        if (det < tolerance)
            return false;

        if (p_point)
            *p_point = (cross(n2, n1) * w) + (cross(n0, n2) * p1.w) - (n0xn1 * p2.w) / det;

        return true;
    }

    GEM_INLINE bool GEM_VECTORCALL plane4f::intersects_ray(const ray3f& ray, float3* p_point, const float tolerance /*= 0.001f*/)
    {
        float3 n = normal();
        float ndotv = dot(n, ray.v);
        if (ndotv < tolerance)
            return false;

        if (p_point)
            *p_point = ray.p - ray.v * (dot(n, ray.p) / ndotv);

        return true;
    }
}