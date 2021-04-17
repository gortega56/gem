#pragma once
#include "common/defines.h"
#include "quaternion.h"

namespace gem
{
#pragma region transform3f

    struct transform3f
    {
        quatf rotation;
        float3 translation;
        float3 scale;

        static transform3f identity();

        transform3f(const quatf& q, const float3& t, const float3& s);

        transform3f(const transform3f& o);

        transform3f() = default;

        transform3f& GEM_VECTORCALL concatenate(const transform3f rhs);

        float4x3 matrix4x3();

        float4x4 matrix4x4();

        transform3f& GEM_VECTORCALL operator=(const transform3f& o);
    };

    transform3f GEM_VECTORCALL transform_mul3f(const transform3f& lhs, const transform3f& rhs);

    GEM_INLINE transform3f transform3f::identity()
    {
        return transform3f(quatf::identity(), float3(0, 0, 0), float3(1, 1, 1));
    }

    GEM_INLINE transform3f::transform3f(const quatf& q, const float3& t, const float3& s)
        : rotation(q)
        , translation(t)
        , scale(s)
    {

    }

    GEM_INLINE transform3f::transform3f(const transform3f& o)
        : rotation(o.rotation)
        , translation(o.translation)
        , scale(o.scale)
    {

    }

    GEM_INLINE transform3f& GEM_VECTORCALL transform3f::concatenate(const transform3f rhs)
    {
        quatf q = rotation * rhs.rotation;
        float3 t = rhs.translation + (rhs.rotation * (translation * rhs.scale));
        float3 s = scale * rhs.scale;

        rotation = q;
        translation = t;
        scale = s;
        
        return *this;
    }

    GEM_INLINE float4x3 transform3f::matrix4x3()
    {
        return rotation.matrix4x3() * float4x3::scale(scale) * float4x3::translate(translation);
    }

    GEM_INLINE float4x4 transform3f::matrix4x4()
    {
        return rotation.matrix4x4() * float4x4::scale(scale) * float4x4::translate(translation);
    }

    GEM_INLINE transform3f& GEM_VECTORCALL transform3f::operator=(const transform3f& o)
    {
        rotation = o.rotation;
        translation = o.translation;
        scale = o.scale;
        return *this;
    }

    GEM_INLINE transform3f GEM_VECTORCALL transform_mul3f(const transform3f& lhs, const transform3f& rhs)
    {
        return
        {
            lhs.rotation * rhs.rotation,
            rhs.translation + (rhs.rotation * (lhs.translation * rhs.scale)),
            lhs.scale * rhs.scale
        };
    }

#pragma endregion

#pragma region transform1f

    struct transform1f
    {
        quatf rotation;
        float3 translation;
        float scale;

        static transform1f identity();

        transform1f(const quatf& q, const float3& t, const float s);

        transform1f(const transform1f& o);

        transform1f() = default;

        transform1f& GEM_VECTORCALL concatenate(const transform1f rhs);

        float4x3 matrix4x3();

        float4x4 matrix4x4();

        transform1f& GEM_VECTORCALL operator=(const transform1f & o);
    };

    transform1f GEM_VECTORCALL transform_mul1f(const transform1f& lhs, const transform1f& rhs);

    GEM_INLINE transform1f transform1f::identity()
    {
        return transform1f(quatf::identity(), float3(0, 0, 0), 1);
    }

    GEM_INLINE transform1f::transform1f(const quatf& q, const float3& t, const float s)
        : rotation(q)
        , translation(t)
        , scale(s)
    {

    }

    GEM_INLINE transform1f::transform1f(const transform1f& o)
        : rotation(o.rotation)
        , translation(o.translation)
        , scale(o.scale)
    {

    }

    GEM_INLINE transform1f& GEM_VECTORCALL transform1f::concatenate(const transform1f rhs)
    {
        quatf q = rotation * rhs.rotation;
        float3 t = rhs.translation + (rhs.rotation * (translation * rhs.scale));
        float s = scale * rhs.scale;

        rotation = q;
        translation = t;
        scale = s;

        return *this;
    }

    GEM_INLINE float4x3 transform1f::matrix4x3()
    {
        return rotation.matrix4x3() * float4x3::scale(scale) * float4x3::translate(translation);
    }

    GEM_INLINE float4x4 transform1f::matrix4x4()
    {
        return rotation.matrix4x4() * float4x4::scale(scale) * float4x4::translate(translation);
    }

    GEM_INLINE transform1f& GEM_VECTORCALL transform1f::operator=(const transform1f& o)
    {
        rotation = o.rotation;
        translation = o.translation;
        scale = o.scale;
        return *this;
    }

    GEM_INLINE transform1f GEM_VECTORCALL transform_mul1f(const transform1f& lhs, const transform1f& rhs)
    {
        return
        {
            lhs.rotation * rhs.rotation,
            rhs.translation + (rhs.rotation * (lhs.translation * rhs.scale)),
            lhs.scale * rhs.scale
        };
    }

#pragma endregion
}