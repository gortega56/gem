#pragma once
#include "line.h"
#include "ray.h"

#include <cfloat>
#include <cmath>

namespace gem
{
    struct plane4f
    {
        float3 n;
        float  d;

        static plane4f implicit(float a, float b, float c, float d);

        static plane4f GEM_VECTORCALL transform(const transform3f& transform, const plane4f& plane);

        static plane4f GEM_VECTORCALL transform(const transform1f& transform, const plane4f& plane);

        static plane4f GEM_VECTORCALL transform(const float4x4& m, const plane4f& plane);

        static plane4f GEM_VECTORCALL transform(const float4x3& m, const plane4f& plane);

        static plane4f GEM_VECTORCALL transform(const float3x3& m, const plane4f& plane);

        void normalize();

        float distance(const float3& point) const;
        
        bool GEM_VECTORCALL intersects(const plane4f& plane, line3f* out, float tolerance = 0.001f) const;

        bool GEM_VECTORCALL intersects(const plane4f& p0, const plane4f& p1, float3* out, float tolerance = 0.001f) const;

        bool GEM_VECTORCALL intersects(const ray3f& ray, float3* out, const float tolerance = 0.001f);
    };

    GEM_INLINE float GEM_VECTORCALL dot(const plane4f& lhs, const float3& rhs);

    GEM_INLINE plane4f GEM_VECTORCALL normalize(const plane4f& p);

    GEM_INLINE plane4f plane4f::implicit(float a, float b, float c, float d)
    {
        float l = gem::length(a, b, c);
        if (l < FLT_EPSILON)
        {
            a = 0;
            b = 0;
            c = 0;
            d = 0;
        }
        else
        {
            float i = (1.f / l);
            a = a * i;
            b = b * i;
            c = c * i;
            d = d * i;
        }
       
        return { { a, b, c }, -d };
    }

    GEM_INLINE plane4f GEM_VECTORCALL plane4f::transform(const transform3f& transform, const plane4f& plane)
    {
        const float3& s = transform.s;
        float3 n = transform.q.transform_point(plane.n * float3{ 1.f / s.x, 1.f / s.y, 1.f / s.z });
        float l = n.length();
        float d = plane.d + dot(n, transform.t);
        return { n / l, d / l };
    }

    GEM_INLINE plane4f GEM_VECTORCALL plane4f::transform(const transform1f& transform, const plane4f& plane)
    {
        float3 n = transform.transform_vector(plane.n);
        float l = n.length();
        float d = plane.d + dot(n, transform.t);
        return { n / l, d / l};
    }

    GEM_INLINE plane4f GEM_VECTORCALL plane4f::transform(const float4x4& m, const plane4f& plane)
    {
        // TODO(gortega): Replace with adj
        float4x4 adj = m.inverse().transpose();
        float4 h = float4(plane.n, 0) * adj;
        float3 n = { h.x, h.y, h.z };
        float3 t = { m.m30, m.m31, m.m32 };
        float l = n.length();
        float d = plane.d + dot(t, n);
        return { n / l, d / l };
        
    }

    GEM_INLINE plane4f GEM_VECTORCALL plane4f::transform(const float4x3& m, const plane4f& plane)
    {
        float4x3 adj = m.inverse().transpose();
        float3 n = float4(plane.n, 0) * adj;
        float3 t = { m.m30, m.m31, m.m32 };
        float l = n.length();
        float d = plane.d + dot(t, n);
        return { n / l, d / l };
    }

    GEM_INLINE plane4f GEM_VECTORCALL plane4f::transform(const float3x3& m, const plane4f& plane)
    {
        float3x3 adj = m.inverse().transpose();
        float3 n = plane.n * adj;
        float l = n.length();
        float d = plane.d;
        return { n / l, d / l };
    }

    GEM_INLINE void plane4f::normalize()
    {
        float l = length(n);
        float il = 1.f / l;
        n *= il;
        d *= il;
    }

    GEM_INLINE float plane4f::distance(const float3& point) const
    {
        return dot(*this, point);
    }

    GEM_INLINE bool GEM_VECTORCALL plane4f::intersects(const plane4f& plane, line3f* p_line, float tolerance /*= 0.001f*/) const
    {
        float3 n0 = n;
        float3 n1 = plane.n;
        float3 v = cross(n0, n1);
        float det = length_squared(v);
        if (det < tolerance)
            return false;

        // Plucker moment of the intersection line, m = q x v, in this file's n.p - d = 0
        // convention. Lengyel's [n | d] form (n.p + d = 0) yields the negation of this.
        if (p_line)
            *p_line = { v, (plane.d * n0) - (d * n1) };

        return true;
    }

    GEM_INLINE bool GEM_VECTORCALL plane4f::intersects(const plane4f& p1, const plane4f& p2, float3* p_point, float tolerance /*= 0.001f*/) const
    {
        float3 n0 = n;
        float3 n1 = p1.n;
        float3 n2 = p2.n;
        float3 n0xn1 = cross(n0, n1);
        float det = dot(n0xn1, n2);
        if (fabsf(det) < tolerance)
            return false;

        // Cramer's rule on [n0; n1; n2] p = [d0; d1; d2].
        if (p_point)
            *p_point = ((cross(n1, n2) * d) + (cross(n2, n0) * p1.d) + (n0xn1 * p2.d)) / det;

        return true;
    }

    GEM_INLINE bool GEM_VECTORCALL plane4f::intersects(const ray3f& ray, float3* p_point, const float tolerance /*= 0.001f*/)
    {
        float fv = dot(n.x, n.y, n.z, ray.v.x, ray.v.y, ray.v.z);
        if (abs(fv) < tolerance)
            return false;

        float fp = dot(*this, ray.p);
        if (p_point)
            *p_point = ray.p - ray.v * (fp / fv);

        return true;
    }

    GEM_INLINE plane4f GEM_VECTORCALL normalize(const plane4f& p)
    {
        float l = length(p.n);
        float il = 1.f / l;
        return { p.n * il, p.d * il };
    }

    GEM_INLINE float GEM_VECTORCALL dot(const plane4f& lhs, const float3& rhs)
    {
        return (lhs.n.x * rhs.x) + (lhs.n.y * rhs.y) + (lhs.n.z * rhs.z) - lhs.d;
    }
}