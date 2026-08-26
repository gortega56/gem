#pragma once
#include "common/defines.h"
#include "quaternion.h"

namespace gem
{
#pragma region affine3f

    /// <summary>
    /// Affine transform of the form: M + t, where t is a translation vector and M is a linear transform.    
    /// Transforms can be concatenated in the form (A + a)B + b = AB + (aB + b).
    /// Points are transformed via pM + t.
    /// </summary>
    struct affine3f
    {
        quatf  q = quatf::identity();
        float3 t = float3::zero();
        float3 s = float3::one();

        static affine3f identity() { return { }; };

        affine3f inverse() const;

        void GEM_VECTORCALL concat(const affine3f& b);

        float3 GEM_VECTORCALL mulv(const float3 p) const;

        float3 GEM_VECTORCALL mulp(const float3 p) const;

        float4x3 matrix4x3() const;

        float4x4 matrix4x4() const;

        affine3f& GEM_VECTORCALL operator=(const affine3f& o);
    };

    affine3f GEM_VECTORCALL concat(const affine3f& a, const affine3f& b);

    GEM_INLINE affine3f GEM_VECTORCALL concat(const affine3f& a, const affine3f& b)
    {
        return {
            quatf::concat(a.q, b.q),
            b.q.mul(a.t * b.s) + b.t,
            a.s * b.s
        };
    }

    GEM_INLINE affine3f affine3f::inverse() const
    {
        quatf invq = gem::inverse(q);
        float3 invs = { 1.f / s.x, 1.f / s.y, 1.f / s.z };
        float3 invt = -(invq.mul(t) * invs);
        return { invq, invt, invs };
    }

    GEM_INLINE void GEM_VECTORCALL affine3f::concat(const affine3f& b)
    {
       *this = gem::concat(*this, b);
    }

    GEM_INLINE float3 GEM_VECTORCALL affine3f::mulv(const float3 p) const
    {
        return q.mul(p * s);
    }

    GEM_INLINE float3 GEM_VECTORCALL affine3f::mulp(const float3 p) const
    {
        return q.mul(p * s) + t;
    }

    GEM_INLINE float4x3 affine3f::matrix4x3() const
    {
        return q.matrix4x3() * float4x3::scale(s) * float4x3::translate(t);
    }

    GEM_INLINE float4x4 affine3f::matrix4x4() const
    {
        return q.matrix4x4() * float4x4::scale(s) * float4x4::translate(t);
    }

    GEM_INLINE affine3f& GEM_VECTORCALL affine3f::operator=(const affine3f& o)
    {
        q = o.q;
        t = o.t;
        s = o.s;
        return *this;
    }

#pragma endregion

#pragma region similarity3f

    /// <summary>
    /// Affine transform of the form: M + t, where t is a translation vector and M is a linear transform with rotation and uniform scale.    
    /// Transforms can be concatenated in the form (A + a)B + b = AB + (aB + b).
    /// Points are transformed via pM + t.
    /// </summary>
    struct similarity3f
    {
        quatf  q = quatf::identity();
        float3 t = float3::zero();
        float  s = 1.f;

        static similarity3f identity() { return { }; };

        similarity3f inverse() const;

        void GEM_VECTORCALL concat(const similarity3f& b);

        float3 GEM_VECTORCALL mulv(const float3 p) const;

        float3 GEM_VECTORCALL mulp(const float3 p) const;

        float4x3 matrix4x3() const;

        float4x4 matrix4x4() const;

        similarity3f& GEM_VECTORCALL operator=(const similarity3f& o);

        affine3f affine() const { return { q, t, { s, s, s } }; }
    };

    similarity3f GEM_VECTORCALL concat(const similarity3f& a, const similarity3f& b);

    GEM_INLINE similarity3f GEM_VECTORCALL concat(const similarity3f& a, const similarity3f& b)
    {
        return {
            quatf::concat(a.q, b.q),
            b.q.mul(a.t * b.s) + b.t,
            a.s * b.s
        };
    }

    GEM_INLINE similarity3f similarity3f::inverse() const
    {
        quatf invq = gem::inverse(q);
        float invs = 1.f / s;
        float3 invt = -(invq.mul(t) * invs);
        return { invq, invt, invs };
    }

    GEM_INLINE void GEM_VECTORCALL similarity3f::concat(const similarity3f& b)
    {
        *this = gem::concat(*this, b);
    }

    GEM_INLINE float3 GEM_VECTORCALL similarity3f::mulv(const float3 p) const
    {
        return q.mul(p * s);
    }

    GEM_INLINE float3 GEM_VECTORCALL similarity3f::mulp(const float3 p) const
    {
        return q.mul(p * s) + t;
    }

    GEM_INLINE float4x3 similarity3f::matrix4x3() const
    {
        return q.matrix4x3() * float4x3::scale(s) * float4x3::translate(t);
    }

    GEM_INLINE float4x4 similarity3f::matrix4x4() const
    {
        return q.matrix4x4() * float4x4::scale(s) * float4x4::translate(t);
    }

    GEM_INLINE similarity3f& GEM_VECTORCALL similarity3f::operator=(const similarity3f& o)
    {
        q = o.q;
        t = o.t;
        s = o.s;
        return *this;
    }

#pragma endregion

#pragma region rigid3f

    /// <summary>
    /// Affine transform of the form: M + t, where t is a translation vector and M is a linear transform with rotation and implicit uniform scale of 1.    
    /// Transforms can be concatenated in the form (A + a)B + b = AB + (aB + b).
    /// Points are transformed via pM + t.
    /// </summary>
    struct rigid3f
    {
        quatf  q = quatf::identity();
        float3 t = float3::zero();

        rigid3f inverse() const;

        void GEM_VECTORCALL concat(const rigid3f& b);

        float3 GEM_VECTORCALL mulv(const float3 p) const;

        float3 GEM_VECTORCALL mulp(const float3 p) const;

        float4x3 matrix4x3() const;

        float4x4 matrix4x4() const;

        rigid3f& GEM_VECTORCALL operator=(const rigid3f& o);
    };

    rigid3f GEM_VECTORCALL concat(const rigid3f& a, const rigid3f& b);

    GEM_INLINE rigid3f GEM_VECTORCALL concat(const rigid3f& a, const rigid3f& b)
    {
        return {
            quatf::concat(a.q, b.q),
            b.q.mul(a.t) + b.t,
        };
    }

    GEM_INLINE rigid3f rigid3f::inverse() const
    {
        quatf invq = gem::inverse(q);
        float3 invt = -invq.mul(t);
        return { invq, invt };
    }

    GEM_INLINE void GEM_VECTORCALL rigid3f::concat(const rigid3f& b)
    {
        *this = gem::concat(*this, b);
    }

    GEM_INLINE float3 GEM_VECTORCALL rigid3f::mulv(const float3 p) const
    {
        return q.mul(p);
    }

    GEM_INLINE float3 GEM_VECTORCALL rigid3f::mulp(const float3 p) const
    {
        return q.mul(p) + t;
    }

    GEM_INLINE float4x3 rigid3f::matrix4x3() const
    {
        return q.matrix4x3() * float4x3::translate(t);
    }

    GEM_INLINE float4x4 rigid3f::matrix4x4() const
    {
        return q.matrix4x4() * float4x4::translate(t);
    }

    GEM_INLINE rigid3f& GEM_VECTORCALL rigid3f::operator=(const rigid3f& o)
    {
        q = o.q;
        t = o.t;
        return *this;
    }

#pragma endregion
}