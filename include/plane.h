#pragma once
#include "line.h"

namespace gem
{
    struct plane4f
    {
        union
        {
            struct
            {
                float x, y, z, w;
            };

            float3 n;
            float d;
        };

        plane4f& normalize();

        float GEM_VECTORCALL distance_from_origin() const;

        bool GEM_VECTORCALL intersects_plane(const plane4f& plane, line3f* p_line, float tolerance = 0.001f) const;

        bool GEM_VECTORCALL intersects_planes(const plane4f& p0, const plane4f& p1, float3* p_point, float tolerance = 0.001f) const;
    };

    GEM_INLINE plane4f& plane4f::normalize()
    {
        float il = 1.0f / sqrtf((x * x) + (y * y) + (z * z));
        x *= il;
        y *= il;
        z *= il;
        w *= il;
        return *this;
    }

    GEM_INLINE float plane4f::distance_from_origin() const
    {
        return fabsf(d) / length(n);
    }

    GEM_INLINE bool GEM_VECTORCALL plane4f::intersects_plane(const plane4f& plane, line3f* p_line, float tolerance /*= 0.001f*/) const
    {
        float3 v = cross(n, plane.n);
        float det = length_squared(v);
        if (det < tolerance)
            return false;

        if (p_line)
            *p_line = { v, (d * plane.n) - (plane.d - n) };

        return true;
    }

    GEM_INLINE bool GEM_VECTORCALL plane4f::intersects_planes(const plane4f& p1, const plane4f& p2, float3* p_point, float tolerance /*= 0.001f*/) const
    {
        float3 n0xn1 = cross(n, p1.n);
        float det = dot(n0xn1, p2.n);
        if (det < tolerance)
            return false;

        if (p_point)
            *p_point = (cross(p2.n, p1.n) * d) + (cross(n, p2.n) * p1.d) - (n0xn1 * p2.d) / det;

        return true;
    }
}