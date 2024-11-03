#pragma once
#include "common/defines.h"
#include "quaternion.h"

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

        transform3f();

        transform3f(const quatf& q, const float3& t, const float3& s);

        transform3f inverse() const;

        transform3f& GEM_VECTORCALL concatenate(const transform3f& b);

        float3 GEM_VECTORCALL transform_vector(const float3& p) const;

        float3 GEM_VECTORCALL transform_point(const float3& p) const;

        float4x3 matrix4x3() const;

        float4x4 matrix4x4() const;

        transform3f& GEM_VECTORCALL operator=(const transform3f& o);
    };

    transform3f GEM_VECTORCALL concatenate(const transform3f& a, const transform3f& b);

    GEM_INLINE transform3f transform3f::identity()
    {
        return transform3f::set(
            quatf::identity(), 
            float3(0, 0, 0), 
            float3(1, 1, 1) 
        );
    }

    GEM_INLINE transform3f GEM_VECTORCALL transform3f::set(const quatf& q, const float3& t, const float3& s)
    {
        return {
            q,
            t,
            s
        };
    }

    GEM_INLINE transform3f GEM_VECTORCALL transform3f::set(const transform3f& o)
    {
        return {
            o.q,
            o.t,
            o.s
        };
    }

    GEM_INLINE transform3f::transform3f()
        : q(quatf::identity())
        , t({0.f, 0.f, 0.f})
        , s({1.f, 1.f, 1.f})
    {

    }

    GEM_INLINE transform3f::transform3f(const quatf& q, const float3& t, const float3& s)
        : q(q)
        , t(t)
        , s(s)
    {

    }

    GEM_INLINE transform3f transform3f::inverse() const
    {
        quatf invq = gem::inverse(q);
        float3 invs = { 1.f / s.x, 1.f / s.y, 1.f / s.z };
        float3 invt = -(invq.transform_point(t) * invs);
        return transform3f::set(invq, invt, invs);
    }

    GEM_INLINE transform3f& GEM_VECTORCALL transform3f::concatenate(const transform3f& b)
    {
        // (A + a)B + b = AB + (aB + b)
        t = b.q.transform_point(t * b.s) + b.t;   // (aB + b)
        q = q * b.q;
        s = s * b.s;
        
        return *this;
    }

    GEM_INLINE float3 GEM_VECTORCALL transform3f::transform_vector(const float3& p) const
    {
        return q.transform_point(p * s);
    }

    GEM_INLINE float3 GEM_VECTORCALL transform3f::transform_point(const float3& p) const
    {
        return q.transform_point(p * s) + t;
    }

    GEM_INLINE float4x3 transform3f::matrix4x3() const
    {
        return q.matrix4x3() * float4x3::scale(s) * float4x3::translate(t);
    }

    GEM_INLINE float4x4 transform3f::matrix4x4() const
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
        return transform3f::set(
            a.q * b.q,
            b.q.transform_point(a.t * b.s) + b.t,
            a.s * b.s
        );
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

        transform1f inverse() const;

        transform1f& GEM_VECTORCALL concatenate(const transform1f& b);

        float3 GEM_VECTORCALL transform_vector(const float3& p) const;

        float3 GEM_VECTORCALL transform_point(const float3& p) const;

        float4x3 matrix4x3() const;

        float4x4 matrix4x4() const;

        transform1f& GEM_VECTORCALL operator=(const transform1f & o);
    };

    transform1f GEM_VECTORCALL concatenate(const transform1f& a, const transform1f& b);

    GEM_INLINE transform1f transform1f::identity()
    {
        return transform1f::set(
            quatf::identity(),
            float3(0, 0, 0),
            1
        );
    }

    GEM_INLINE transform1f GEM_VECTORCALL transform1f::set(const quatf& q, const float3& t, const float s)
    {
        return transform1f::set(
            q,
            t,
            s
        );
    }

    GEM_INLINE transform1f GEM_VECTORCALL transform1f::set(const transform1f& o)
    {
        return transform1f::set(
            o.q,
            o.t,
            o.s
        );
    }

    GEM_INLINE transform1f transform1f::inverse() const
    {
        quatf invq = gem::inverse(q);
        float invs = { 1.f / s };
        float3 invt = -(invq.transform_point(t) * invs);
        return  transform1f::set(invq, invt, invs );
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

    GEM_INLINE float3 GEM_VECTORCALL transform1f::transform_vector(const float3& p) const
    {
        return q.transform_point(p) * s;
    }

    GEM_INLINE float3 GEM_VECTORCALL transform1f::transform_point(const float3& p) const
    {
        return (q.transform_point(p) * s) + t;
    }

    GEM_INLINE float4x3 transform1f::matrix4x3() const
    {
        return q.matrix4x3() * float4x3::scale(s) * float4x3::translate(t);
    }

    GEM_INLINE float4x4 transform1f::matrix4x4() const
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
        return transform1f::set(
            a.q * b.q,
            b.t + (b.q.transform_point(a.t * b.s)),
            a.s * b.s
        );
    }

#pragma endregion
}