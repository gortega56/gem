#pragma once
#include "common/defines.h"
#include "quaternion.h"
#include "ray.h"
#include "line.h"
#include "plane.h"

namespace gem
{
#pragma region transform3f  

    /// <summary>
    /// Affine transform of the form: M + t, where t is a translation vector and M is a linear transform.    
    /// Transforms can be concatenated in the form (A + a)B + b = AB + (aB + b).
    /// Points are transformed via pM + t.
    /// </summary>
    struct transform3f
    {
        quatf q;
        float3 t;
        float3 s;

        static transform3f identity();

        static transform3f GEM_VECTORCALL set(const quatf& q, const float3& t, const float3& s);

        static transform3f GEM_VECTORCALL set(const transform3f& o);

        transform3f() = default;

        transform3f& GEM_VECTORCALL concatenate(const transform3f& b);

        float3 GEM_VECTORCALL transform_vector(const float3& p);

        float3 GEM_VECTORCALL transform_point(const float3& p);

        ray3f GEM_VECTORCALL transform_ray(const ray3f& p);

        line3f GEM_VECTORCALL transform_line(const line3f& p);

        plane4f GEM_VECTORCALL transform_plane(const plane4f& p);

        float4x3 matrix4x3();

        float4x4 matrix4x4();

        transform3f& GEM_VECTORCALL operator=(const transform3f& o);
    };

    transform3f GEM_VECTORCALL concatenate(const transform3f& a, const transform3f& b);

    GEM_INLINE transform3f transform3f::identity()
    {
        return
        {
            quatf::identity(), 
            float3(0, 0, 0), 
            float3(1, 1, 1) 
        };
    }

    GEM_INLINE transform3f GEM_VECTORCALL transform3f::set(const quatf& q, const float3& t, const float3& s)
    {
        return
        {
            q,
            t,
            s
        };
    }

    GEM_INLINE transform3f GEM_VECTORCALL transform3f::set(const transform3f& o)
    {
        return
        {
            o.q,
            o.t,
            o.s
        };
    }

    GEM_INLINE transform3f& GEM_VECTORCALL transform3f::concatenate(const transform3f& b)
    {
        quatf  _q = q * b.q;
        float3 _t = b.t + (b.q.transform_point(t * b.s));
        float3 _s = s * b.s;

        q = _q;
        t = _t;
        s = _s;
        
        return *this;
    }

    GEM_INLINE float3 GEM_VECTORCALL transform3f::transform_vector(const float3& p)
    {
        return q.transform_point(p) * s;
    }

    GEM_INLINE float3 GEM_VECTORCALL transform3f::transform_point(const float3& p)
    {
        return (q.transform_point(p) * s) + t;
    }

    GEM_INLINE ray3f GEM_VECTORCALL transform3f::transform_ray(const ray3f& p)
    {
        return { transform_point(p.p), transform_vector(p.v) };
    }

    GEM_INLINE line3f GEM_VECTORCALL transform3f::transform_line(const line3f& p)
    {
        float4x3 h = matrix4x3().inverse().transpose();
        return { transform_vector(p.v), p.m * h };
    }

    GEM_INLINE plane4f GEM_VECTORCALL transform3f::transform_plane(const plane4f& p)
    {
        float4x4 m = matrix4x4().inverse().transpose();
        float4 q = { p.x, p.y, p.z, p.d };
        float4 r = q * m;
        return { r.x, r.y, r.z, r.w };
    }

    GEM_INLINE float4x3 transform3f::matrix4x3()
    {
        return q.matrix4x3() * float4x3::scale(s) * float4x3::translate(t);
    }

    GEM_INLINE float4x4 transform3f::matrix4x4()
    {
        return q.matrix4x4() * float4x4::scale(s) * float4x4::translate(t);
    }

    GEM_INLINE transform3f& GEM_VECTORCALL transform3f::operator=(const transform3f& o)
    {
        q = o.q;
        t = o.t;
        s = o.s;
        return *this;
    }

    GEM_INLINE transform3f GEM_VECTORCALL concatenate(const transform3f& a, const transform3f& b)
    {
        return
        {
            a.q * b.q,
            b.t + (b.q.transform_point(a.t * b.s)),
            a.s * b.s
        };
    }

#pragma endregion

#pragma region transform1f

    /// <summary>
    /// Affine transform of the form: M + t, where t is a translation vector and M is a linear transform with rotation and uniform scale.    
    /// Transforms can be concatenated in the form (A + a)B + b = AB + (aB + b).
    /// Points are transformed via pM + t.
    /// </summary>
    struct transform1f
    {
        quatf q;
        float3 t;
        float s;

        static transform1f identity();

        static transform1f GEM_VECTORCALL set(const quatf& q, const float3& t, const float s);

        static transform1f GEM_VECTORCALL set(const transform1f& o);

        transform1f() = default;

        transform1f& GEM_VECTORCALL concatenate(const transform1f& b);

        float3 GEM_VECTORCALL transform_vector(const float3& p);

        float3 GEM_VECTORCALL transform_point(const float3& p);

        ray3f GEM_VECTORCALL transform_ray(const ray3f& p);

        line3f GEM_VECTORCALL transform_line(const line3f& p);

        plane4f GEM_VECTORCALL transform_plane(const plane4f& p);

        float4x3 matrix4x3();

        float4x4 matrix4x4();

        transform1f& GEM_VECTORCALL operator=(const transform1f & o);
    };

    transform1f GEM_VECTORCALL concatenate(const transform1f& a, const transform1f& b);

    GEM_INLINE transform1f transform1f::identity()
    {
        return
        {
            quatf::identity(),
            float3(0, 0, 0),
            1
        };
    }

    GEM_INLINE transform1f GEM_VECTORCALL transform1f::set(const quatf& q, const float3& t, const float s)
    {
        return
        {
            q,
            t,
            s
        };
    }

    GEM_INLINE transform1f GEM_VECTORCALL transform1f::set(const transform1f& o)
    {
        return
        {
            o.q,
            o.t,
            o.s
        };
    }

    GEM_INLINE transform1f& GEM_VECTORCALL transform1f::concatenate(const transform1f& b)
    {
        quatf  _q = q * b.q;
        float3 _t = b.t + (b.q.transform_point(t * b.s));
        float  _s = s * b.s;

        q = _q;
        t = _t;
        s = _s;

        return *this;
    }

    GEM_INLINE float3 GEM_VECTORCALL transform1f::transform_vector(const float3& p)
    {
        return q.transform_point(p) * s;
    }

    GEM_INLINE float3 GEM_VECTORCALL transform1f::transform_point(const float3& p)
    {
        return (q.transform_point(p) * s) + t;
    }

    GEM_INLINE ray3f GEM_VECTORCALL transform1f::transform_ray(const ray3f& p)
    {
        return { transform_point(p.p), transform_vector(p.v) };
    }

    GEM_INLINE line3f GEM_VECTORCALL transform1f::transform_line(const line3f& p)
    {
        float4x3 h = matrix4x3().inverse().transpose();
        return { transform_vector(p.v), p.m * h };
    }

    GEM_INLINE plane4f GEM_VECTORCALL transform1f::transform_plane(const plane4f& p)
    {
        float4x4 m = matrix4x4().inverse().transpose();
        float4 q = { p.x, p.y, p.z, p.d };
        float4 r = q * m;
        return { r.x, r.y, r.z, r.w };
    }

    GEM_INLINE float4x3 transform1f::matrix4x3()
    {
        return q.matrix4x3() * float4x3::scale(s) * float4x3::translate(t);
    }

    GEM_INLINE float4x4 transform1f::matrix4x4()
    {
        return q.matrix4x4() * float4x4::scale(s) * float4x4::translate(t);
    }

    GEM_INLINE transform1f& GEM_VECTORCALL transform1f::operator=(const transform1f& o)
    {
        q = o.q;
        t = o.t;
        s = o.s;
        return *this;
    }

    GEM_INLINE transform1f GEM_VECTORCALL concatenate(const transform1f& a, const transform1f& b)
    {
        return
        {
            a.q * b.q,
            b.t + (b.q.transform_point(a.t * b.s)),
            a.s * b.s
        };
    }

#pragma endregion
}