#pragma once
#include "transform.h"

#define PLUCKER 0

namespace gem
{
    struct line2f
    {
        float a, b, c;

        float2 GEM_VECTORCALL point_closest_to_point(const float2& point) const;

        float GEM_VECTORCALL distance_to_point(const float2& point) const;

        bool GEM_VECTORCALL intersects_line(const line2f& line, float2* p_point, const float tolerance = 0.001f) const;
    };

    GEM_INLINE float2 GEM_VECTORCALL line2f::point_closest_to_point(const float2& p) const
    {
        float2 n = { a, b };
        float il = 1.0f / length(n);
        n *= il;
        float2 o = n * (c * il);
        float d = length(o - project(p, o));
        return p + (n * d);
    }

    GEM_INLINE float GEM_VECTORCALL line2f::distance_to_point(const float2& p) const
    {
        float2 n = { a, b };
        float il = 1.0f / length(n);
        n *= il;
        float2 o = n * (c * il);
        return length(o - project(p, o));
    }

    GEM_INLINE bool GEM_VECTORCALL line2f::intersects_line(const line2f& line, float2* p_point, const float tolerance /*= 0.001f*/) const
    {
        float det = (a * line.b) - (line.a * b);
        if (det < tolerance)
            return false;

        if (p_point)
            *p_point = { (c * line.b) - (line.c * b), (a * line.c) - (c * line.a) };

        return true;
    }

    #pragma region plucker
   
    struct plucker3f
    {
        float3 v; // direction
        float3 m; // moment
    
        static float GEM_VECTORCALL distance(const float3& v0, const float3& m0, const float3& v1, const float3& m1, const float tolerance = 0.001f);

        static float GEM_VECTORCALL distance(const float3& v, const float3& m, const float3& p);

        static bool GEM_VECTORCALL intersects(const float3& v0, const float3& m0, const float3& v1, const float3& m1, float3* p, const float threshold, const float tolerance = 0.001f);

        static plucker3f GEM_VECTORCALL transform(const plucker3f& l, const affine3f& transform);

        static plucker3f GEM_VECTORCALL transform(const plucker3f& l, const similarity3f& transform);

        static plucker3f GEM_VECTORCALL transform(const plucker3f& l, const float4x4& m);

        static plucker3f GEM_VECTORCALL transform(const plucker3f& l, const float4x3& m);

        static plucker3f GEM_VECTORCALL transform(const plucker3f& l, const float3x3& m);

        static plucker3f GEM_VECTORCALL make(const float3 p0, const float3 p1);
    
        float GEM_VECTORCALL distance(const plucker3f& l, const float tolerance = 0.001f);

        float GEM_VECTORCALL distance(const float3& p);

        bool GEM_VECTORCALL intersects(const plucker3f& l, float3* p, const float threshold, const float tolerance = 0.001f);

        void GEM_VECTORCALL transform(const affine3f& transform);

        void GEM_VECTORCALL transform(const similarity3f& transform);

        void GEM_VECTORCALL transform(const float4x4& m);

        void GEM_VECTORCALL transform(const float4x3& m);

        void GEM_VECTORCALL transform(const float3x3& m);
    };

    GEM_INLINE plucker3f GEM_VECTORCALL plucker3f::make(const float3 p0, const float3 p1)
    {
        return { p1 - p0, cross(p0, p1) };
    }

    GEM_INLINE float GEM_VECTORCALL plucker3f::distance(const float3& v0, const float3& m0, const float3& v1, const float3& m1, const float tolerance /*= 0.001f*/)
    {
        float l0 = length(v0);
        float l1 = length(v1);
        float v0xv1 = length(cross(v0, v1));
        if (v0xv1 <= (l0 * l1 * tolerance))
        {
            float3 q = cross(v1, m1) / dot(v1, v1);
            return distance(v0, m0, q);
        }

        return fabsf(dot(v0, m1) + dot(v1, m0)) / v0xv1;
    }

    GEM_INLINE float GEM_VECTORCALL plucker3f::distance(const float3& v, const float3& m, const float3& p)
    {
        return length(cross(v, p) + m) / length(v);
    }

    GEM_INLINE bool GEM_VECTORCALL plucker3f::intersects(const float3& v0, const float3& m0, const float3& v1, const float3& m1, float3* p, const float threshold, const float tolerance /*= 0.001f*/)
    {
        if (distance(v0, m0, v1, m1, tolerance) < threshold)
        {
            if (p)
            {
                float3 v0xv1 = cross(v0, v1);
                float nn = dot(v0xv1, v0xv1);
                float lim = length(v0) * length(v1) * tolerance;
                if (nn <= (lim * lim))
                {
                    return false;
                }

                float3 a = cross(v0, m0) / dot(v0, v0);
                float3 b = cross(v1, m1) / dot(v1, v1);
                *p = a + v0 * (dot(cross(b - a, v1), v0xv1) / nn);
            }

            return true;
        }

        return false;
    }

    GEM_INLINE plucker3f GEM_VECTORCALL plucker3f::transform(const plucker3f& l, const affine3f& x)
    {
        return transform(l, x.matrix4x3());
    }

    GEM_INLINE plucker3f GEM_VECTORCALL plucker3f::transform(const plucker3f& l, const similarity3f& x)
    {
        return transform(l, x.matrix4x3());
    }

    GEM_INLINE plucker3f GEM_VECTORCALL plucker3f::transform(const plucker3f& l, const float4x4& m)
    {
        float4 vM = float4(l.v, 0) * m;
        float4 mA = float4(l.m, 0) * m.adj();
        float3 t = { m.m30, m.m31, m.m32 };
        float3 v = { vM.x, vM.y, vM.z };
        float3 n = { mA.x, mA.y, mA.z };
        return {
            v, 
            n + cross(t, v)
        };
    }

    GEM_INLINE plucker3f GEM_VECTORCALL plucker3f::transform(const plucker3f& l, const float4x3& m)
    {
        float3 vM = l.v * m;
        float3 t = { m.m30, m.m31, m.m32 };
        return {
            vM,
            l.m * m.adj() + cross(t, vM)
        };
    }

    GEM_INLINE plucker3f GEM_VECTORCALL plucker3f::transform(const plucker3f& l, const float3x3& m)
    {
        return {
            l.v * m,
            l.m * m.adj()
        };
    }

    GEM_INLINE float GEM_VECTORCALL plucker3f::distance(const plucker3f& l, const float tolerance /*= 0.001f*/)
    {
        return distance(v, m, l.v, l.m, tolerance);
    }

    GEM_INLINE float GEM_VECTORCALL plucker3f::distance(const float3& p)
    {
        return distance(v, m, p);
    }

    GEM_INLINE bool GEM_VECTORCALL plucker3f::intersects(const plucker3f& l, float3* p, const float threshold, const float tolerance /*= 0.001f*/)
    {
        return intersects(v, m, l.v, l.m, p, threshold, tolerance);
    }

    GEM_INLINE void GEM_VECTORCALL plucker3f::transform(const affine3f& x)
    {
        *this = transform(*this, x);
    }

    GEM_INLINE void GEM_VECTORCALL plucker3f::transform(const similarity3f& x)
    {
        *this = transform(*this, x);
    }

    GEM_INLINE void GEM_VECTORCALL plucker3f::transform(const float4x4& m)
    {
        *this = transform(*this, m);
    }

    GEM_INLINE void GEM_VECTORCALL plucker3f::transform(const float4x3& m)
    {
        *this = transform(*this, m);
    }

    GEM_INLINE void GEM_VECTORCALL plucker3f::transform(const float3x3& m)
    {
        *this = transform(*this, m);
    }

    #pragma endregion

    #pragma region line3f

    struct line3f
    {
        float3 p0, p1;

        static line3f GEM_VECTORCALL closest(const float3& p0, const float3& p1, const float3& q0, const float3& q1, const float tolerance = 0.001f);

        static float2 GEM_VECTORCALL closest_t(const float3& p0, const float3& p1, const float3& q0, const float3& q1, const float tolerance = 0.001f);

        static float GEM_VECTORCALL closest_t(const float3& p0, const float3& p1, const float3& q);

        static float GEM_VECTORCALL distance(const float3& p0, const float3& p1, const float3& q0, const float3& q1, const float tolerance = 0.001f);

        static float GEM_VECTORCALL distance(const float3& p0, const float3& p1, const float3& q);

        static bool GEM_VECTORCALL intersects(const float3& p0, const float3& p1, const float3& q0, const float3& q1, float3* p, const float threshold, const float tolerance = 0.001f);

        static line3f GEM_VECTORCALL transform(const line3f& l, const affine3f& transform);

        static line3f GEM_VECTORCALL transform(const line3f& l, const similarity3f& transform);

        static line3f GEM_VECTORCALL transform(const line3f& l, const float4x4& m);

        static line3f GEM_VECTORCALL transform(const line3f& l, const float4x3& m);

        static line3f GEM_VECTORCALL transform(const line3f& l, const float3x3& m);

        line3f GEM_VECTORCALL closest(const line3f& l, const float tolerance = 0.001f);

        float3 GEM_VECTORCALL closest(const float3& p);

        float GEM_VECTORCALL distance(const line3f& l, const float tolerance = 0.001f);

        float GEM_VECTORCALL distance(const float3& p);

        bool GEM_VECTORCALL intersects(const line3f& l, float3* p, const float threshold, const float tolerance = 0.001f);

        void GEM_VECTORCALL transform(const affine3f& transform);

        void GEM_VECTORCALL transform(const similarity3f& transform);

        void GEM_VECTORCALL transform(const float4x4& m);

        void GEM_VECTORCALL transform(const float4x3& m);

        void GEM_VECTORCALL transform(const float3x3& m);
    };

    GEM_INLINE line3f GEM_VECTORCALL line3f::closest(const float3& p0, const float3& p1, const float3& q0, const float3& q1, const float tolerance /*=0.001f*/)
    {
        float2 t = closest_t(p0, p1, q0, q1, tolerance);
        return 
        {
           gem::lerp(p0, p1, t.x),
           gem::lerp(q0, q1, t.y)
        };
    }

    GEM_INLINE float2 GEM_VECTORCALL line3f::closest_t(const float3& p0, const float3& p1, const float3& q0, const float3& q1, const float tolerance /*=0.001f*/)
    {
        float3 d1 = p1 - p0;
        float3 d2 = q1 - q0;
        float3 r = p0 - q0;
        float a = dot(d1,d1);
        float b = dot(d1,d2);
        float e = dot(d2,d2);
        float d = (a*e) - (b*b);
        if (d <= (a * e * tolerance))   // sin^2(theta) <= tolerance : parallel (scale invariant)
        {
            float s = 0;
            float t = (e > 0.f) ? closest_t(q0, q1, p0) : 0.f;
            return { s, t };
        }

        float c = dot(d1, r);
        float f = dot(d2, r);
        float s = ((b*f)-(c*e))/d;
        float t = ((a*f)-(b*c))/d;
        return { s, t };
        /*float3 v0 = p1 - p0;
        float3 v1 = q1 - q0;
        float3 w = q0 - p0;
        float m0 = v0.length_squared();
        float m1 = v1.length_squared();
        float v0v1 = dot(v0, v1);
        float s0 = dot(w, v0);
        float s1 = dot(w, v1);
        float d = (v0v1 * v0v1) - (m0 * m1);
        if (d > -tolerance)
        {
            return {
                0.0f,
                -s1 / m1
            };
        }

        d = 1.f / d;
        float t0 = d * ((v0v1 * s1) - (m1 * s0));
        float t1 = d * ((m0 * s1) - (v0v1 * s0));

        return { t0, t1 };*/
    }

    GEM_INLINE float GEM_VECTORCALL line3f::closest_t(const float3& p0, const float3& p1, const float3& q)
    {
        float3 ab = p1 - p0;
        float3 ac = q  - p0;
        return dot(ac, ab) / dot(ab, ab);
    }

    GEM_INLINE float GEM_VECTORCALL line3f::distance(const float3& p0, const float3& p1, const float3& q)
    {
    #if PLUCKER
        return plucker3f::distance(p1 - p0, cross(p0, p1), q);
    #else
        float t = closest_t(p0, p1, q);
        return length(q - lerp(p0, p1, t));
    #endif
    }

    GEM_INLINE float GEM_VECTORCALL line3f::distance(const float3& p0, const float3& p1, const float3& q0, const float3& q1, const float tolerance /*=0.001f*/)
    {
#if PLUCKER
        return plucker3f::distance(p1 - p0, cross(p0, p1), q1 - q0, cross(q0, q1), tolerance);
#else
        float3 d1 = p1 - p0;
        float3 d2 = q1 - q0;
        float3 d1xd2 = cross(d1, d2);
        float l1 = length(d1);
        float l2 = length(d2);
        float denom = length(d1xd2);
        if (denom <= (l1 * l2 * tolerance))   // sin(theta) <= tolerance : parallel (scale invariant)
        {
            return length(cross(q0 - p0, d1)) / l1;
        }

        return fabsf(dot(q0 - p0, d1xd2)) / denom;
#endif
    }

    GEM_INLINE bool GEM_VECTORCALL line3f::intersects(const float3& p0, const float3& p1, const float3& q0, const float3& q1, float3* p, const float threshold, const float tolerance /*=0.001f*/)
    {
#if PLUCKER
        return plucker3f::intersects(p1 - p0, cross(p0, p1), q1 - q0, cross(q0, q1), p, threshold, tolerance);
#else
        line3f cp = closest(p0, p1, q0, q1, tolerance);
        if (length_squared(cp.p1 - cp.p0) < threshold * threshold)
        {
            *p = (cp.p0 + cp.p1) * 0.5f;
            return true;
        }

        return false;
#endif
    }

    GEM_INLINE line3f GEM_VECTORCALL line3f::transform(const line3f& l, const affine3f& transform)
    {
        return {
            transform.mulp(l.p0),
            transform.mulp(l.p1)
        };
    }

    GEM_INLINE line3f GEM_VECTORCALL line3f::transform(const line3f& l, const similarity3f& transform)
    {
        return {
            transform.mulp(l.p0),
            transform.mulp(l.p1)
        };
    }

    GEM_INLINE line3f GEM_VECTORCALL line3f::transform(const line3f& l, const float4x4& m)
    {
        float4 p0 = float4(l.p0, 1) * m;
        float4 p1 = float4(l.p1, 1) * m;
        return {
            { p0.x, p0.y, p0.z },
            { p1.x, p1.y, p1.z }
        };
    }

    GEM_INLINE line3f GEM_VECTORCALL line3f::transform(const line3f& l, const float4x3& m)
    {
        return {
            l.p0 * m,
            l.p1 * m
        };
    }
    
    GEM_INLINE line3f GEM_VECTORCALL line3f::transform(const line3f& l, const float3x3& m)
    {
        return {
            l.p0 * m,
            l.p1 * m
        };
    }

    GEM_INLINE line3f GEM_VECTORCALL line3f::closest(const line3f& l, const float tolerance /*=0.001f*/)
    {
        return closest(p0, p1, l.p0, l.p1, tolerance);
    }

    GEM_INLINE float3 GEM_VECTORCALL line3f::closest(const float3& p)
    {
        float3 ab = p1 - p0;
        float3 ac = p  - p0;
        float t = dot(ac, ab) / dot(ab, ab);
        return p0 + ab * t;
    }

    GEM_INLINE float GEM_VECTORCALL line3f::distance(const line3f& l, const float tolerance /*=0.001f*/)
    {
        return distance(p0, p1, l.p0, l.p1, tolerance);
    }

    GEM_INLINE float GEM_VECTORCALL line3f::distance(const float3& p)
    {
        return distance(p0, p1, p);
    }

    GEM_INLINE bool GEM_VECTORCALL line3f::intersects(const line3f& l, float3* p, const float threshold, const float tolerance /*=0.001f*/)
    {
        return intersects(p0, p1, l.p0, l.p1, p, threshold, tolerance);
    }

    GEM_INLINE void GEM_VECTORCALL line3f::transform(const affine3f& x)
    {
        *this = transform(*this, x.matrix4x3());
    }

    GEM_INLINE void GEM_VECTORCALL line3f::transform(const similarity3f& x)
    {
        *this = transform(*this, x.matrix4x3());
    }

    GEM_INLINE void GEM_VECTORCALL line3f::transform(const float4x4& m)
    {
        *this = transform(*this, m);
    }

    GEM_INLINE void GEM_VECTORCALL line3f::transform(const float4x3& m)
    {
        *this = transform(*this, m);
    }

    GEM_INLINE void GEM_VECTORCALL line3f::transform(const float3x3& m)
    {
       *this = transform(*this, m);
    }

    #pragma endregion

    struct segment2f
    {
        float2 a, b;

        float2 GEM_VECTORCALL point_closest_to_point(const float2& point) const;
        
        float GEM_VECTORCALL distance_to_point(const float2& point) const;

        bool GEM_VECTORCALL intersects_line(const segment2f& line, float2* p_point, const float tolerance = 0.001f) const;
    };

    GEM_INLINE float2 GEM_VECTORCALL segment2f::point_closest_to_point(const float2& point) const
    {
        return project(point - a, b - a);
    }

    GEM_INLINE float GEM_VECTORCALL segment2f::distance_to_point(const float2& point) const
    {
        float d = (point.x * a.y) - (point.y * a.x);
        float l = length_squared(b - a);
        return sqrtf((d * d) / l);
    }

    GEM_INLINE bool GEM_VECTORCALL segment2f::intersects_line(const segment2f& line, float2* p_point, const float tolerance /*= 0.001f*/) const
    {
        float2 v0 = b - a;
        float2 v1 = line.b - line.a;
        float2 dp = line.a - a;
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
        if (t < 0)
            return false;

        if (t > 1)
            return false;

        if (p_point)
            *p_point = line.a + v0 * t;
        
        return true;
    }

    struct segment3f
    {
        float3 p0, p1;

        static segment3f GEM_VECTORCALL closest(const float3 p0, const float3 p1, const float3 q0, const float3 q1, const float tolerance = 0.001f);

        static float3 GEM_VECTORCALL closest(const float3 p0, const float3 p1, const float3 p);

        static float2 GEM_VECTORCALL closest_t(const float3 p0, const float3 p1, const float3 q0, const float3 q1, const float tolerance = 0.001f);
    
        static float GEM_VECTORCALL closest_t(const float3 p0, const float3 p1, const float3 p);

        static float GEM_VECTORCALL distance(const float3 p0, const float3 p1, const float3 q0, const float3 q1, const float tolerance = 0.001f);

        static float GEM_VECTORCALL distance(const float3 p0, const float3 p1, const float3 p);

        static bool GEM_VECTORCALL intersects(const float3 p0, const float3 p1, const float3 q0, const float3 q1, float3* p, const float threshold, const float tolerance = 0.001f);

        static segment3f GEM_VECTORCALL transform(const segment3f& l, const affine3f& transform);

        static segment3f GEM_VECTORCALL transform(const segment3f& l, const similarity3f& transform);

        static segment3f GEM_VECTORCALL transform(const segment3f& l, const float4x4& m);

        static segment3f GEM_VECTORCALL transform(const segment3f& l, const float4x3& m);

        static segment3f GEM_VECTORCALL transform(const segment3f& l, const float3x3& m);

        segment3f GEM_VECTORCALL closest(const segment3f& l, const float tolerance = 0.001f);

        float3 GEM_VECTORCALL closest(const float3 p);

        float GEM_VECTORCALL distance(const segment3f& l, const float tolerance = 0.001f);

        float GEM_VECTORCALL distance(const float3 p);

        bool GEM_VECTORCALL intersects(const float3 q0, const float3 q1, float3* p, const float threshold, const float tolerance = 0.001f);

        void GEM_VECTORCALL transform(const affine3f& m);

        void GEM_VECTORCALL transform(const similarity3f& m);

        void GEM_VECTORCALL transform(const float4x4& m);

        void GEM_VECTORCALL transform(const float4x3& m);

        void GEM_VECTORCALL transform(const float3x3& m);
    };

    GEM_INLINE segment3f GEM_VECTORCALL segment3f::closest(const float3 p0, const float3 p1, const float3 q0, const float3 q1, const float tolerance /*=0.001f*/)
    {
        float2 t = closest_t(p0, p1, q0, q1, tolerance);
        return {
            lerp(p0, p1, t.x),
            lerp(q0, q1, t.y)
        };
    }

    GEM_INLINE float3 GEM_VECTORCALL segment3f::closest(const float3 p0, const float3 p1, const float3 p)
    {
        float t = closest_t(p0, p1, p);
        return lerp(p0, p1, t);
    }

    GEM_INLINE float clamp01(float v)
    {
        if (v < 0.f) return 0.f;
        if (v > 1.f) return 1.f;
        return v;
    }

    GEM_INLINE float2 GEM_VECTORCALL segment3f::closest_t(const float3 p0, const float3 p1, const float3 q0, const float3 q1, const float tolerance /*=0.001f*/)
    {
        float3 d1 = p1 - p0;
        float3 d2 = q1 - q0;
        float3 r = p0 - q0;
        float a = dot(d1,d1);
        float e = dot(d2,d2);
        float f = dot(d2, r);
        float tol = tolerance * tolerance;
        if (a <= tol && e <= tol)
        {
            return { 0, 0 };
        }

        float s, t;
        if (a <= tol)
        {
            s = 0.f;
            t = f / e;
            t = clamp01(t);
        }
        else
        {
            float c = dot(d1, r);
            if (e <= tol)
            {
                t = 0.f;
                s = clamp01(-c / a);
            }
            else
            {
                float b = dot(d1,d2);
                float denom = (a*e) - (b*b);
                if (denom != 0.f)
                {
                    s = clamp01(((b*f) - (c*e)) / denom);
                }
                else
                {
                    s = 0.f;
                }

                t = ((b*s) + f) / e;

                if (t < 0.f)
                {
                    t = 0.f;
                    s = clamp01(-c/a);
                }
                else if (t > 1.f)
                {
                    t = 1.f;
                    s = clamp01((b-c)/a);
                }
            }
        }

        return { s, t };
    }

    GEM_INLINE float GEM_VECTORCALL segment3f::closest_t(const float3 p0, const float3 p1, const float3 p)
    {
        float3 ab = p1 - p0;
        float3 ac = p  - p0;
        return clamp01(dot(ac, ab)/dot(ab,ab));
    }

    GEM_INLINE float GEM_VECTORCALL segment3f::distance(const float3 p0, const float3 p1, const float3 q0, const float3 q1, const float tolerance /*=0.001f*/)
    {
        segment3f s = closest(p0, p1, q0, q1, tolerance);
        return length(s.p1 - s.p0);
    }

    GEM_INLINE float GEM_VECTORCALL segment3f::distance(const float3 p0, const float3 p1, const float3 p)
    {
        float t = closest_t(p0, p1, p);
        return length(p - lerp(p0, p1, t));
    }

    GEM_INLINE bool GEM_VECTORCALL segment3f::intersects(const float3 p0, const float3 p1, const float3 q0, const float3 q1, float3* p, const float threshold, const float tolerance /*=0.001f*/)
    {
        segment3f cp = closest(p0, p1, q0, q1, tolerance);
        if (length_squared(cp.p1 - cp.p0) < threshold * threshold)
        {
            *p = (cp.p0 + cp.p1) * 0.5f;
            return true;
        }

        return false;
    }

    GEM_INLINE segment3f GEM_VECTORCALL segment3f::transform(const segment3f& l, const affine3f& m)
    {
        return {
                m.mulp(l.p0),
                m.mulp(l.p1)
        };
    }

    GEM_INLINE segment3f GEM_VECTORCALL segment3f::transform(const segment3f& l, const similarity3f& m)
    {
        return {
            m.mulp(l.p0),
            m.mulp(l.p1)
        };
    }

    GEM_INLINE segment3f GEM_VECTORCALL segment3f::transform(const segment3f& l, const float4x4& m)
    {
        float4 p0 = float4(l.p0, 1) * m;
        float4 p1 = float4(l.p1, 1) * m;
        return {
            { p0.x, p0.y, p0.z },
            { p1.x, p1.y, p1.z }
        };
    }

    GEM_INLINE segment3f GEM_VECTORCALL segment3f::transform(const segment3f& l, const float4x3& m)
    {
        return { l.p0 * m, l.p1 * m };
    }

    GEM_INLINE segment3f GEM_VECTORCALL segment3f::transform(const segment3f& l, const float3x3& m)
    {
        return { l.p0 * m, l.p1 * m };
    }

    GEM_INLINE segment3f GEM_VECTORCALL segment3f::closest(const segment3f& l, const float tolerance /*=0.001f*/)
    {
        return closest(p0, p1, l.p0, l.p1, tolerance);
    }

    GEM_INLINE float3 GEM_VECTORCALL segment3f::closest(const float3 p)
    {
        return closest(p0, p1, p);
    }

    GEM_INLINE float GEM_VECTORCALL segment3f::distance(const segment3f& l, const float tolerance /*=0.001f*/)
    {
        return distance(p0, p1, l.p0, l.p1, tolerance);
    }

    GEM_INLINE float GEM_VECTORCALL segment3f::distance(const float3 p)
    {
        return distance(p0, p1, p);
    }

    GEM_INLINE bool GEM_VECTORCALL segment3f::intersects(const float3 q0, const float3 q1, float3* p, const float threshold, const float tolerance /*=0.001f*/)
    {
        return intersects(p0, p1, q0, q1, p, threshold, tolerance);
    }

    GEM_INLINE void GEM_VECTORCALL segment3f::transform(const affine3f& m)
    {
        *this = transform(*this, m);
    }

    GEM_INLINE void GEM_VECTORCALL segment3f::transform(const similarity3f& m)
    {
        *this = transform(*this, m);
    }

    GEM_INLINE void GEM_VECTORCALL segment3f::transform(const float4x4& m)
    {
        *this = transform(*this, m);
    }

    GEM_INLINE void GEM_VECTORCALL segment3f::transform(const float4x3& m)
    {
        *this = transform(*this, m);
    }

    GEM_INLINE void GEM_VECTORCALL segment3f::transform(const float3x3& m)
    {
        *this = transform(*this, m);
    }
}