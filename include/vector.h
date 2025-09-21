#pragma once
#include "common/defines.h"
#include <cmath>

namespace gem
{
#pragma region float2

    struct float2
    {
        float x, y;

        float2(const float ix, const float iy);

        float2(const float* p_components);

        float2(const float component);

        float2(const float2& o);

        float2() = default;

        float length_squared() const;

        float length() const;

        float2& normalize();

        float2& safe_normalize(const float tolerance = 0.001f);

        float2& GEM_VECTORCALL operator+=(const float2& rhs);

        float2& GEM_VECTORCALL operator-=(const float2& rhs);

        float2& GEM_VECTORCALL operator*=(const float2& rhs);

        float2& GEM_VECTORCALL operator=(const float2& rhs);

        float2& operator*=(const float rhs);

        float2& operator/=(const float rhs);

        float& operator[](const unsigned int index);

        float operator[](const unsigned int index) const;

        float2 operator-() const;
    };

    float GEM_VECTORCALL dot(const float lx, const float ly, const float rx, const float ry);

    float GEM_VECTORCALL dot(const float2& lhs, const float2& rhs);

    float GEM_VECTORCALL length_squared(const float x, const float y);

    float GEM_VECTORCALL length_squared(const float2& val);

    float GEM_VECTORCALL length(const float x, const float y);

    float GEM_VECTORCALL length(const float2& rhs);

    float2 GEM_VECTORCALL normalize(const float2& rhs);

    float2 GEM_VECTORCALL safe_normalize(const float2& rhs, const float tolerance = 0.001f);

    float2 GEM_VECTORCALL lerp(const float2& lhs, const float2& rhs, const float t);

    float2 GEM_VECTORCALL project(const float2& lhs, const float2& rhs);

    float2 GEM_VECTORCALL reject(const float2& lhs, const float2& rhs);

    float2 GEM_VECTORCALL operator+(const float2& lhs, const float2& rhs);

    float2 GEM_VECTORCALL operator-(const float2& lhs, const float2& rhs);

    float2 GEM_VECTORCALL operator*(const float2& lhs, const float2& rhs);

    float2 GEM_VECTORCALL operator*(const float2& lhs, const float rhs);

    float2 GEM_VECTORCALL operator/(const float2& lhs, const float rhs);

    float2 GEM_VECTORCALL operator*(const float lhs, const float2& rhs);

    GEM_INLINE float2::float2(const float ix, const float iy)
        : x(ix)
        , y(iy)
    {

    }

    GEM_INLINE float2::float2(const float* p_components)
        : x(p_components[0])
        , y(p_components[1])
    {

    }

    GEM_INLINE float2::float2(const float component)
        : x(component)
        , y(component)
    {

    }

    GEM_INLINE float2::float2(const float2& o)
        : x(o.x)
        , y(o.y)
    {

    }

    GEM_INLINE float float2::length_squared() const
    {
        return (x * x) + (y * y);
    }

    GEM_INLINE float float2::length() const
    {
        return sqrtf((x * x) + (y * y));
    }

    GEM_INLINE float2& float2::normalize()
    {
        float il = 1.0f / sqrtf((x * x) + (y * y));
        x *= il;
        y *= il;
        return *this;
    }

    GEM_INLINE float2& float2::safe_normalize(const float tolerance /*= 0.001f*/)
    {
        float l = sqrtf((x * x) + (y * y));
        if (l > tolerance)
        {
            float il = 1.0f / l;
            x *= il;
            y *= il;
        }

        return *this;
    }

    GEM_INLINE float2& GEM_VECTORCALL float2::operator+=(const float2& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    GEM_INLINE float2& GEM_VECTORCALL float2::operator-=(const float2& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    GEM_INLINE float2& GEM_VECTORCALL float2::operator*=(const float2& rhs)
    {
        x *= rhs.x;
        y *= rhs.y;
        return *this;
    }

    GEM_INLINE float2& GEM_VECTORCALL float2::operator=(const float2& rhs)
    {
        x = rhs.x;
        y = rhs.y;
        return *this;
    }

    GEM_INLINE float2& float2::operator*=(const float rhs)
    {
        x *= rhs;
        y *= rhs;
        return *this;
    }

    GEM_INLINE float2& float2::operator/=(const float rhs)
    {
        float inv = 1.0f / rhs;
        x *= inv;
        y *= inv;
        return *this;
    }

    GEM_INLINE float& float2::operator[](const unsigned int index)
    {
        return reinterpret_cast<float*>(this)[index];
    }

    GEM_INLINE float float2::operator[](const unsigned int index) const
    {
        return reinterpret_cast<const float*>(this)[index];
    }

    GEM_INLINE float2 float2::operator-() const
    {
        return
        {
            -x,
            -y
        };
    }

    GEM_INLINE float GEM_VECTORCALL dot(const float lx, const float ly, const float rx, const float ry)
    {
        return (lx * rx) + (ly * ry);
    }

    GEM_INLINE float GEM_VECTORCALL dot(const float2& lhs, const float2& rhs)
    {
        return (lhs.x * rhs.x) + (lhs.y * rhs.y);
    }

    GEM_INLINE float GEM_VECTORCALL length_squared(const float x, const float y)
    {
        return sqrtf((x * x) + (y * y));
    }

    GEM_INLINE float GEM_VECTORCALL length_squared(const float2& val)
    {
        return (val.x * val.x) + (val.y * val.y);
    }

    GEM_INLINE float GEM_VECTORCALL length(const float x, const float y)
    {
        return sqrtf((x * x) + (y * y));
    }

    GEM_INLINE float GEM_VECTORCALL length(const float2& rhs)
    {
        return sqrtf((rhs.x * rhs.x) + (rhs.y * rhs.y));
    }

    GEM_INLINE float2 GEM_VECTORCALL normalize(const float2& rhs)
    {
        float il = 1.0f / sqrtf((rhs.x * rhs.x) + (rhs.y * rhs.y));
        return
        {
            rhs.x * il,
            rhs.y * il
        };
    }

    GEM_INLINE float2 GEM_VECTORCALL safe_normalize(const float2& rhs, const float tolerance /*= 0.001f*/)
    {
        float l = sqrtf((rhs.x * rhs.x) + (rhs.y * rhs.y));
        float il = (l > tolerance) ? (1.0f / l) : 1.0f;
        return
        {
            rhs.x * il,
            rhs.y * il
        };
    }

    GEM_INLINE float2 GEM_VECTORCALL lerp(const float2& lhs, const float2& rhs, const float t)
    {
        float oneMinusT = 1.0f - t;
        return
        {
            (oneMinusT * lhs.x) + (rhs.x * t),
            (oneMinusT * lhs.y) + (rhs.y * t)
        };
    }

    GEM_INLINE float2 GEM_VECTORCALL project(const float2& lhs, const float2& rhs)
    {
        float d = (lhs.x * rhs.x) + (lhs.y * rhs.y);
        float l = (rhs.x * rhs.x) + (rhs.y * rhs.y);
        float t = d / l;
        return
        {
            rhs.x * t,
            rhs.y * t
        };
    }

    GEM_INLINE float2 GEM_VECTORCALL reject(const float2& lhs, const float2& rhs)
    {
        float d = (lhs.x * rhs.x) + (lhs.y * rhs.y);
        float l = (rhs.x * rhs.x) + (rhs.y * rhs.y);
        float t = d / l;
        return
        {
            lhs.x - rhs.x * t,
            lhs.y - rhs.y * t
        };
    }

    GEM_INLINE float2 GEM_VECTORCALL operator+(const float2& lhs, const float2& rhs)
    {
        return
        {
            lhs.x + rhs.x,
            lhs.y + rhs.y
        };
    }

    GEM_INLINE float2 GEM_VECTORCALL operator-(const float2& lhs, const float2& rhs)
    {
        return
        {
            lhs.x - rhs.x,
            lhs.y - rhs.y
        };
    }

    GEM_INLINE float2 GEM_VECTORCALL operator*(const float2& lhs, const float2& rhs)
    {
        return
        {
            lhs.x * rhs.x,
            lhs.y * rhs.y
        };
    }

    GEM_INLINE float2 GEM_VECTORCALL operator*(const float2& lhs, const float rhs)
    {
        return
        {
            lhs.x * rhs,
            lhs.y * rhs
        };
    }

    GEM_INLINE float2 GEM_VECTORCALL operator/(const float2& lhs, const float rhs)
    {
        float inv = 1.0f / rhs;
        return
        {
            lhs.x * inv,
            lhs.y * inv
        };
    }

    GEM_INLINE float2 GEM_VECTORCALL operator*(const float lhs, const float2& rhs)
    {
        return
        {
            rhs.x * lhs,
            rhs.y * lhs
        };
    }

#pragma endregion

#pragma region float3

    struct float3
    {
        float x, y, z;

        static float3 zero();

        static float3 one();

        static float3 init(const float x, const float y, const float z);

        static float3 init(const float xyz);

        static float3 init(const float* xyz);

        static float3 init(const float2 xy, const float z);

        float length_squared() const;

        float length() const;

        float3& normalize();

        float3& safe_normalize(const float tolerance = 0.001f);

        float3& GEM_VECTORCALL operator+=(const float3& rhs);

        float3& GEM_VECTORCALL operator-=(const float3& rhs);

        float3& GEM_VECTORCALL operator*=(const float3& rhs);

        float3& operator*=(const float rhs);

        float3& operator/=(const float rhs);

        float3 operator-() const;

        float& operator[](const unsigned int index);

        float operator[](const unsigned int index) const;
    };

    float GEM_VECTORCALL dot(const float lx, const float ly, const float lz, const float rx, const float ry, const float rz);

    float GEM_VECTORCALL dot(const float3& lhs, const float3& rhs);

    float GEM_VECTORCALL length_squared(const float x, const float y, const float z);

    float GEM_VECTORCALL length_squared(const float3& val);

    float GEM_VECTORCALL length(const float x, const float y, const float z);

    float GEM_VECTORCALL length(const float3& rhs);

    float3 GEM_VECTORCALL normalize(const float3& rhs);

    float3 GEM_VECTORCALL safe_normalize(const float3& rhs, const float tolerance = 0.001f);

    float3 GEM_VECTORCALL lerp(const float3& lhs, const float3& rhs, const float t);

    float3 GEM_VECTORCALL project(const float3& lhs, const float3& rhs);

    float3 GEM_VECTORCALL reject(const float3& lhs, const float3& rhs);

    float3 GEM_VECTORCALL rotate_axis_angle(const float3& v, const float3& a, const float angle);

    float3 GEM_VECTORCALL cross(const float3& lhs, const float3& rhs);

    float3 GEM_VECTORCALL operator+(const float3& lhs, const float3& rhs);

    float3 GEM_VECTORCALL operator-(const float3& lhs, const float3& rhs);

    float3 GEM_VECTORCALL operator*(const float3& lhs, const float3& rhs);

    float3 GEM_VECTORCALL operator*(const float3& lhs, const float rhs);

    float3 GEM_VECTORCALL operator/(const float3& lhs, const float rhs);

    float3 GEM_VECTORCALL operator*(const float lhs, const float3& rhs);

    GEM_INLINE float3 float3::zero()
    {
        return { 0, 0, 0 };
    }

    GEM_INLINE float3 float3::one()
    {
        return { 1, 1, 1 };
    }

    GEM_INLINE float3 float3::init(const float x, const float y, const float z)
    {
        return { x, y, z };
    }

    GEM_INLINE float3 float3::init(const float xyz)
    {
        return { xyz, xyz, xyz };
    }

    GEM_INLINE float3 float3::init(const float* xyz)
    {
        return { xyz[0], xyz[1], xyz[2] };
    }

    GEM_INLINE float3 float3::init(const float2 xy, const float z)
    {
        return { xy.x, xy.y, z };
    }

    GEM_INLINE float float3::length_squared() const
    {
        return (x * x) + (y * y) + (z * z);
    }

    GEM_INLINE float float3::length() const
    {
        return sqrtf((x * x) + (y * y) + (z * z));
    }

    GEM_INLINE float3& float3::normalize()
    {
        float il = 1.0f / sqrtf((x * x) + (y * y) + (z * z));
        x *= il;
        y *= il;
        z *= il;
        return *this;
    }

    GEM_INLINE float3& float3::safe_normalize(const float tolerance /*= 0.001f*/)
    {
        float l = sqrtf((x * x) + (y * y) + (z * z));
        if (l > tolerance)
        {
            float il = 1.0f / l;
            x *= il;
            y *= il;
            z *= il;
        }

        return *this;
    }

    GEM_INLINE float3& GEM_VECTORCALL float3::operator+=(const float3& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    GEM_INLINE float3& GEM_VECTORCALL float3::operator-=(const float3& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    GEM_INLINE float3& GEM_VECTORCALL float3::operator*=(const float3& rhs)
    {
        x *= rhs.x;
        y *= rhs.y;
        z *= rhs.z;
        return *this;
    }

    GEM_INLINE float3& float3::operator*=(const float rhs)
    {
        x *= rhs;
        y *= rhs;
        z *= rhs;
        return *this;
    }

    GEM_INLINE float3& float3::operator/=(const float rhs)
    {
        float inv = 1.0f / rhs;
        x *= inv;
        y *= inv;
        z *= inv;
        return *this;
    }

    GEM_INLINE float3 float3::operator-() const
    {
        return
        {
            -x,
            -y,
            -z
        };
    }

    GEM_INLINE float& float3::operator[](const unsigned int index)
    {
        return reinterpret_cast<float*>(this)[index];
    }

    GEM_INLINE float float3::operator[](const unsigned int index) const
    {
        return reinterpret_cast<const float*>(this)[index];
    }

    GEM_INLINE float GEM_VECTORCALL dot(const float lx, const float ly, const float lz, const float rx, const float ry, const float rz)
    {
        return (lx * rx) + (ly * ry) + (lz * rz);
    }

    GEM_INLINE float GEM_VECTORCALL dot(const float3& lhs, const float3& rhs)
    {
        return (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z);
    }

    GEM_INLINE float GEM_VECTORCALL length_squared(const float x, const float y, const float z)
    {
        return (x * x) + (y * y) + (z * z);
    }

    GEM_INLINE float GEM_VECTORCALL length_squared(const float3& val)
    {
        return (val.x * val.x) + (val.y * val.y) + (val.z * val.z);
    }

    GEM_INLINE float GEM_VECTORCALL length(const float x, const float y, const float z)
    {
        return sqrtf((x * x) + (y * y) + (z * z));
    }

    GEM_INLINE float GEM_VECTORCALL length(const float3& rhs)
    {
        return sqrtf((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z));
    }

    GEM_INLINE float3 GEM_VECTORCALL normalize(const float3& rhs)
    {
        float il = 1.0f / sqrtf((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z));
        return
        {
            rhs.x * il,
            rhs.y * il,
            rhs.z * il
        };
    }

    GEM_INLINE float3 GEM_VECTORCALL safe_normalize(const float3& rhs, const float tolerance /*= 0.001f*/)
    {
        float l = sqrtf((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z));
        float il = (l > tolerance) ? (1.0f / l) : 1.0f;
        return
        {
            rhs.x * il,
            rhs.y * il,
            rhs.z * il
        };
    }

    GEM_INLINE float3 GEM_VECTORCALL lerp(const float3& lhs, const float3& rhs, const float t)
    {
        float oneMinusT = 1.0f - t;
        return
        {
            (oneMinusT * lhs.x) + (rhs.x * t),
            (oneMinusT * lhs.y) + (rhs.y * t),
            (oneMinusT * lhs.z) + (rhs.z * t)
        };
    }

    GEM_INLINE float3 GEM_VECTORCALL project(const float3& lhs, const float3& rhs)
    {
        float d = (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z);
        float l = (rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z);
        float t = d / l;
        return
        {
            rhs.x * t,
            rhs.y * t,
            rhs.z * t
        };
    }

    GEM_INLINE float3 GEM_VECTORCALL reject(const float3& lhs, const float3& rhs)
    {
        float d = (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z);
        float l = (rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z);
        float t = d / l;
        return
        {
            lhs.x - rhs.x * t,
            lhs.y - rhs.y * t,
            lhs.z - rhs.z * t
        };
    }

    GEM_INLINE float3 GEM_VECTORCALL rotate_axis_angle(const float3& v, const float3& a, const float theta)
    {
        // v' = (v project a) + (v reject a) * cos(theta) + (a x v) * sin(theta)
        float cos_theta = cosf(theta);
        float c1 = (1.0f - cos_theta) * dot(v, a);
        float sin_theta = sinf(theta);
        return (v * cos_theta) + (a * c1) + cross(a, v) * sin_theta;
    }

    GEM_INLINE float3 GEM_VECTORCALL cross(const float3& lhs, const float3& rhs)
    {
        return
        {
            (lhs.y * rhs.z) - (lhs.z * rhs.y),
            (lhs.z * rhs.x) - (lhs.x * rhs.z),
            (lhs.x * rhs.y) - (lhs.y * rhs.x)
        };
    }

    GEM_INLINE float3 GEM_VECTORCALL operator+(const float3& lhs, const float3& rhs)
    {
        return
        {
            lhs.x + rhs.x,
            lhs.y + rhs.y,
            lhs.z + rhs.z
        };
    }

    GEM_INLINE float3 GEM_VECTORCALL operator-(const float3& lhs, const float3& rhs)
    {
        return
        {
            lhs.x - rhs.x,
            lhs.y - rhs.y,
            lhs.z - rhs.z
        };
    }

    GEM_INLINE float3 GEM_VECTORCALL operator*(const float3& lhs, const float3& rhs)
    {
        return
        {
            lhs.x * rhs.x,
            lhs.y * rhs.y,
            lhs.z * rhs.z
        };
    }

    GEM_INLINE float3 GEM_VECTORCALL operator*(const float3& lhs, const float rhs)
    {
        return
        {
            lhs.x * rhs,
            lhs.y * rhs,
            lhs.z * rhs
        };
    }

    GEM_INLINE float3 GEM_VECTORCALL operator/(const float3& lhs, const float rhs)
    {
        float inv = 1.0f / rhs;
        return
        {
            lhs.x * inv,
            lhs.y * inv,
            lhs.z * inv
        };
    }

    GEM_INLINE float3 GEM_VECTORCALL operator*(const float lhs, const float3& rhs)
    {
        return
        {
            rhs.x * lhs,
            rhs.y * lhs,
            rhs.z * lhs
        };
    }

#pragma endregion

#pragma region float4

    struct float4
    {
        float x, y, z, w;

        float4(const float ix, const float iy, const float iz, const float iw);

        float4(const float* p_components);

        float4(const float component);

        float4(const float3& component, const float w);

        float4(const float2& component, const float z, const float w);

        float4(const float4& o);

        float4() = default;

        float length_squared() const;

        float length() const;

        float4& normalize();

        float4& safe_normalize(const float tolerance = 0.001f);

        float4& GEM_VECTORCALL operator+=(const float4& rhs);

        float4& GEM_VECTORCALL operator-=(const float4& rhs);

        float4& GEM_VECTORCALL operator*=(const float4& rhs);

        float4& GEM_VECTORCALL operator=(const float4& rhs);

        float4& operator*=(const float rhs);

        float4& operator/=(const float rhs);

        float& operator[](const unsigned int index);

        float operator[](const unsigned int index) const;

        float4 operator-() const;
    };

    float GEM_VECTORCALL dot(const float lx, const float ly, const float lz, const float lw, const float rx, const float ry, const float rz, const float rw);

    float GEM_VECTORCALL dot(const float4& lhs, const float4& rhs);

    float GEM_VECTORCALL length_squared(const float x, const float y, const float z, const float w);

    float GEM_VECTORCALL length_squared(const float4& val);

    float GEM_VECTORCALL length(const float x, const float y, const float z, const float w);

    float GEM_VECTORCALL length(const float4& rhs);

    float4 GEM_VECTORCALL normalize(const float4& rhs);

    float4 GEM_VECTORCALL safe_normalize(const float4& rhs, const float tolerance = 0.001f);

    float4 GEM_VECTORCALL lerp(const float4& lhs, const float4& rhs, const float t);

    float4 GEM_VECTORCALL project(const float4& lhs, const float4& rhs);

    float4 GEM_VECTORCALL reject(const float4& lhs, const float4& rhs);

    float4 GEM_VECTORCALL operator+(const float4& lhs, const float4& rhs);

    float4 GEM_VECTORCALL operator-(const float4& lhs, const float4& rhs);

    float4 GEM_VECTORCALL operator*(const float4& lhs, const float4& rhs);

    float4 GEM_VECTORCALL operator*(const float4& lhs, const float rhs);

    float4 GEM_VECTORCALL operator/(const float4& lhs, const float rhs);

    float4 GEM_VECTORCALL operator*(const float lhs, const float4& rhs);

    GEM_INLINE float4::float4(const float ix, const float iy, const float iz, const float iw)
        : x(ix)
        , y(iy)
        , z(iz)
        , w(iw)
    {

    }

    GEM_INLINE float4::float4(const float* p_components)
        : x(p_components[0])
        , y(p_components[1])
        , z(p_components[2])
        , w(p_components[3])
    {

    }

    GEM_INLINE float4::float4(const float component)
        : x(component)
        , y(component)
        , z(component)
        , w(component)
    {

    }

    GEM_INLINE float4::float4(const float4& o)
        : x(o.x)
        , y(o.y)
        , z(o.z)
        , w(o.w)
    {

    }

    GEM_INLINE float4::float4(const float3& component, const float w)
        : x(component.x)
        , y(component.y)
        , z(component.z)
        , w(w)
    {

    }

    GEM_INLINE float4::float4(const float2& component, const float z, const float w)
        : x(component.x)
        , y(component.y)
        , z(z)
        , w(w)
    {

    }

    GEM_INLINE float float4::length_squared() const
    {
        return (x * x) + (y * y) + (z * z) + (w * w);
    }

    GEM_INLINE float float4::length() const
    {
        return sqrtf((x * x) + (y * y) + (z * z) + (w * w));
    }

    GEM_INLINE float4& float4::normalize()
    {
        float il = 1.0f / sqrtf((x * x) + (y * y) + (z * z) + (w * w));
        x *= il;
        y *= il;
        z *= il;
        w *= il;
        return *this;
    }

    GEM_INLINE float4& float4::safe_normalize(const float tolerance /*= 0.001f*/)
    {
        float l = sqrtf((x * x) + (y * y) + (z * z) + (w * w));
        if (l > tolerance)
        {
            float il = 1.0f / l;
            x *= il;
            y *= il;
            z *= il;
            w *= il;
        }

        return *this;
    }

    GEM_INLINE float4& GEM_VECTORCALL float4::operator+=(const float4& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        w += rhs.w;
        return *this;
    }

    GEM_INLINE float4& GEM_VECTORCALL float4::operator-=(const float4& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        w -= rhs.w;
        return *this;
    }

    GEM_INLINE float4& GEM_VECTORCALL float4::operator*=(const float4& rhs)
    {
        x *= rhs.x;
        y *= rhs.y;
        z *= rhs.z;
        w *= rhs.w;
        return *this;
    }

    GEM_INLINE float4& GEM_VECTORCALL float4::operator=(const float4& rhs)
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        w = rhs.w;
        return *this;
    }

    GEM_INLINE float4& float4::operator*=(const float rhs)
    {
        x *= rhs;
        y *= rhs;
        z *= rhs;
        w *= rhs;
        return *this;
    }

    GEM_INLINE float4& float4::operator/=(const float rhs)
    {
        float inv = 1.0f / rhs;
        x *= inv;
        y *= inv;
        z *= inv;
        w *= inv;
        return *this;
    }

    GEM_INLINE float& float4::operator[](const unsigned int index)
    {
        return reinterpret_cast<float*>(this)[index];
    }

    GEM_INLINE float float4::operator[](const unsigned int index) const
    {
        return reinterpret_cast<const float*>(this)[index];
    }

    GEM_INLINE float4 float4::operator-() const
    {
        return
        {
            -x,
            -y,
            -z,
            -w
        };
    }

    GEM_INLINE float GEM_VECTORCALL dot(const float lx, const float ly, const float lz, const float lw, const float rx, const float ry, const float rz, const float rw)
    {
        return (lx * rx) + (ly * ry) + (lz * rz) + (lw * rw);
    }

    GEM_INLINE float GEM_VECTORCALL dot(const float4& lhs, const float4& rhs)
    {
        return (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z) + (lhs.w * rhs.w);
    }

    GEM_INLINE float GEM_VECTORCALL length_squared(const float x, const float y, const float z, const float w)
    {
        return (x * x) + (y * y) + (z * z) + (w * w);
    }

    GEM_INLINE float GEM_VECTORCALL length_squared(const float4& val)
    {
        return (val.x * val.x) + (val.y * val.y) + (val.z * val.z) + (val.w * val.w);
    }

    GEM_INLINE float GEM_VECTORCALL length(const float x, const float y, const float z, const float w)
    {
        return sqrtf((x * x) + (y * y) + (z * z) + (w * w));
    }

    GEM_INLINE float GEM_VECTORCALL length(const float4& rhs)
    {
        return sqrtf((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w));
    }

    GEM_INLINE float4 GEM_VECTORCALL normalize(const float4& rhs)
    {
        float il = 1.0f / sqrtf((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w));
        return
        {
            rhs.x * il,
            rhs.y * il,
            rhs.z * il,
            rhs.w * il
        };
    }

    GEM_INLINE float4 GEM_VECTORCALL safe_normalize(const float4& rhs, const float tolerance /*= 0.001f*/)
    {
        float l = sqrtf((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w));
        float il = (l > tolerance) ? (1.0f / l) : 1.0f;
        return
        {
            rhs.x * il,
            rhs.y * il,
            rhs.z * il,
            rhs.w * il
        };
    }

    GEM_INLINE float4 GEM_VECTORCALL lerp(const float4& lhs, const float4& rhs, const float t)
    {
        float oneMinusT = 1.0f - t;
        return
        {
            (oneMinusT * lhs.x) + (rhs.x * t),
            (oneMinusT * lhs.y) + (rhs.y * t),
            (oneMinusT * lhs.z) + (rhs.z * t),
            (oneMinusT * lhs.w) + (rhs.w * t)
        };
    }

    GEM_INLINE float4 GEM_VECTORCALL project(const float4& lhs, const float4& rhs)
    {
        float d = (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z) + (lhs.w * rhs.w);
        float l = (rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w);
        float t = d / l;
        return
        {
            rhs.x * t,
            rhs.y * t,
            rhs.z * t,
            rhs.w * t
        };
    }

    GEM_INLINE float4 GEM_VECTORCALL reject(const float4& lhs, const float4& rhs)
    {
        float d = (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z) + (lhs.w * rhs.w);
        float l = (rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w);
        float t = d / l;
        return
        {
            lhs.x - rhs.x * t,
            lhs.y - rhs.y * t,
            lhs.z - rhs.z * t,
            lhs.w - rhs.w * t
        };
    }

    GEM_INLINE float4 GEM_VECTORCALL operator+(const float4& lhs, const float4& rhs)
    {
        return
        {
            lhs.x + rhs.x,
            lhs.y + rhs.y,
            lhs.z + rhs.z,
            lhs.w + rhs.w
        };
    }

    GEM_INLINE float4 GEM_VECTORCALL operator-(const float4& lhs, const float4& rhs)
    {
        return
        {
            lhs.x - rhs.x,
            lhs.y - rhs.y,
            lhs.z - rhs.z,
            lhs.w - rhs.w
        };
    }

    GEM_INLINE float4 GEM_VECTORCALL operator*(const float4& lhs, const float4& rhs)
    {
        return
        {
            lhs.x * rhs.x,
            lhs.y * rhs.y,
            lhs.z * rhs.z,
            lhs.w * rhs.w
        };
    }

    GEM_INLINE float4 GEM_VECTORCALL operator*(const float4& lhs, const float rhs)
    {
        return
        {
            lhs.x * rhs,
            lhs.y * rhs,
            lhs.z * rhs,
            lhs.w * rhs
        };
    }

    GEM_INLINE float4 GEM_VECTORCALL operator/(const float4& lhs, const float rhs)
    {
        float inv = 1.0f / rhs;
        return
        {
            lhs.x * inv,
            lhs.y * inv,
            lhs.z * inv,
            lhs.w * inv
        };
    }

    GEM_INLINE float4 GEM_VECTORCALL operator*(const float lhs, const float4& rhs)
    {
        return
        {
            rhs.x * lhs,
            rhs.y * lhs,
            rhs.z * lhs,
            rhs.w * lhs
        };
    }

#pragma endregion 

#pragma region double2

    struct double2
    {
        double x, y;

        double2(const double ix, const double iy);

        double2(const double* p_components);

        double2(const double component);

        double2(const double2& o);

        double2() = default;

        double length_squared() const;

        double length() const;

        double2& normalize();

        double2& safe_normalize(const double tolerance = 0.001f);

        double2& GEM_VECTORCALL operator+=(const double2& rhs);

        double2& GEM_VECTORCALL operator-=(const double2& rhs);

        double2& GEM_VECTORCALL operator*=(const double2& rhs);

        double2& GEM_VECTORCALL operator=(const double2& rhs);

        double2& GEM_VECTORCALL operator*=(const double rhs);

        double2& GEM_VECTORCALL operator/=(const double rhs);

        double& operator[](const unsigned int index);

        double operator[](const unsigned int index) const;

        double2 operator-() const;
    };

    double GEM_VECTORCALL dot(const double lx, const double ly, const double rx, const double ry);

    double GEM_VECTORCALL dot(const double2& lhs, const double2& rhs);

    double GEM_VECTORCALL length_squared(const double x, const double y);

    double GEM_VECTORCALL length_squared(const double2& val);

    double GEM_VECTORCALL length(const double x, const double y);

    double GEM_VECTORCALL length(const double2& rhs);

    double2 GEM_VECTORCALL normalize(const double2& rhs);

    double2 GEM_VECTORCALL safe_normalize(const double2& rhs, const double tolerance = 0.001f);

    double2 GEM_VECTORCALL lerp(const double2& lhs, const double2& rhs, const double t);

    double2 GEM_VECTORCALL project(const double2& lhs, const double2& rhs);

    double2 GEM_VECTORCALL reject(const double2& lhs, const double2& rhs);

    double2 GEM_VECTORCALL operator+(const double2& lhs, const double2& rhs);

    double2 GEM_VECTORCALL operator-(const double2& lhs, const double2& rhs);

    double2 GEM_VECTORCALL operator*(const double2& lhs, const double2& rhs);

    double2 GEM_VECTORCALL operator*(const double2& lhs, const double rhs);

    double2 GEM_VECTORCALL operator/(const double2& lhs, const double rhs);

    double2 GEM_VECTORCALL operator*(const double lhs, const double2& rhs);

    GEM_INLINE double2::double2(const double ix, const double iy)
        : x(ix)
        , y(iy)
    {

    }

    GEM_INLINE double2::double2(const double* p_components)
        : x(p_components[0])
        , y(p_components[1])
    {

    }

    GEM_INLINE double2::double2(const double component)
        : x(component)
        , y(component)
    {

    }

    GEM_INLINE double2::double2(const double2& o)
        : x(o.x)
        , y(o.y)
    {

    }

    GEM_INLINE double double2::length_squared() const
    {
        return (x * x) + (y * y);
    }

    GEM_INLINE double double2::length() const
    {
        return sqrt((x * x) + (y * y));
    }

    GEM_INLINE double2& double2::normalize()
    {
        double il = 1.0f / sqrt((x * x) + (y * y));
        x *= il;
        y *= il;
        return *this;
    }

    GEM_INLINE double2& double2::safe_normalize(const double tolerance /*= 0.001f*/)
    {
        double l = sqrt((x * x) + (y * y));
        if (l > tolerance)
        {
            double il = 1.0f / l;
            x *= il;
            y *= il;
        }

        return *this;
    }

    GEM_INLINE double2& GEM_VECTORCALL double2::operator+=(const double2& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    GEM_INLINE double2& GEM_VECTORCALL double2::operator-=(const double2& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    GEM_INLINE double2& GEM_VECTORCALL double2::operator*=(const double2& rhs)
    {
        x *= rhs.x;
        y *= rhs.y;
        return *this;
    }

    GEM_INLINE double2& GEM_VECTORCALL double2::operator=(const double2& rhs)
    {
        x = rhs.x;
        y = rhs.y;
        return *this;
    }

    GEM_INLINE double2& GEM_VECTORCALL double2::operator*=(const double rhs)
    {
        x *= rhs;
        y *= rhs;
        return *this;
    }

    GEM_INLINE double2& GEM_VECTORCALL double2::operator/=(const double rhs)
    {
        double inv = 1.0f / rhs;
        x *= inv;
        y *= inv;
        return *this;
    }

    GEM_INLINE double& double2::operator[](const unsigned int index)
    {
        return reinterpret_cast<double*>(this)[index];
    }

    GEM_INLINE double double2::operator[](const unsigned int index) const
    {
        return reinterpret_cast<const double*>(this)[index];
    }

    GEM_INLINE double2 double2::operator-() const
    {
        return
        {
            -x,
            -y
        };
    }

    GEM_INLINE double GEM_VECTORCALL dot(const double lx, const double ly, const double rx, const double ry)
    {
        return (lx * rx) + (ly * ry);
    }

    GEM_INLINE double GEM_VECTORCALL dot(const double2& lhs, const double2& rhs)
    {
        return (lhs.x * rhs.x) + (lhs.y * rhs.y);
    }

    GEM_INLINE double GEM_VECTORCALL  length_squared(const double x, const double y)
    {
        return sqrt((x * x) + (y * y));
    }

    GEM_INLINE double GEM_VECTORCALL  length_squared(const double2& val)
    {
        return (val.x * val.x) + (val.y * val.y);
    }

    GEM_INLINE double GEM_VECTORCALL length(const double x, const double y)
    {
        return sqrt((x * x) + (y * y));
    }

    GEM_INLINE double GEM_VECTORCALL length(const double2& rhs)
    {
        return sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y));
    }

    GEM_INLINE double2 GEM_VECTORCALL normalize(const double2& rhs)
    {
        double il = 1.0f / sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y));
        return
        {
            rhs.x * il,
            rhs.y * il
        };
    }

    GEM_INLINE double2 GEM_VECTORCALL safe_normalize(const double2& rhs, const double tolerance /*= 0.001f*/)
    {
        double l = sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y));
        double il = (l > tolerance) ? (1.0f / l) : 1.0f;
        return
        {
            rhs.x * il,
            rhs.y * il
        };
    }

    GEM_INLINE double2 GEM_VECTORCALL lerp(const double2& lhs, const double2& rhs, const double t)
    {
        double oneMinusT = 1.0f - t;
        return
        {
            (oneMinusT * lhs.x) + (rhs.x * t),
            (oneMinusT * lhs.y) + (rhs.y * t)
        };
    }

    GEM_INLINE double2 GEM_VECTORCALL project(const double2& lhs, const double2& rhs)
    {
        double d = (lhs.x * rhs.x) + (lhs.y * rhs.y);
        double l = (rhs.x * rhs.x) + (rhs.y * rhs.y);
        double t = d / l;
        return
        {
            rhs.x * t,
            rhs.y * t
        };
    }

    GEM_INLINE double2 GEM_VECTORCALL reject(const double2& lhs, const double2& rhs)
    {
        double d = (lhs.x * rhs.x) + (lhs.y * rhs.y);
        double l = (rhs.x * rhs.x) + (rhs.y * rhs.y);
        double t = d / l;
        return
        {
            lhs.x - rhs.x * t,
            lhs.y - rhs.y * t
        };
    }

    GEM_INLINE double2 GEM_VECTORCALL operator+(const double2& lhs, const double2& rhs)
    {
        return
        {
            lhs.x + rhs.x,
            lhs.y + rhs.y
        };
    }

    GEM_INLINE double2 GEM_VECTORCALL operator-(const double2& lhs, const double2& rhs)
    {
        return
        {
            lhs.x - rhs.x,
            lhs.y - rhs.y
        };
    }

    GEM_INLINE double2 GEM_VECTORCALL operator*(const double2& lhs, const double2& rhs)
    {
        return
        {
            lhs.x * rhs.x,
            lhs.y * rhs.y
        };
    }

    GEM_INLINE double2 GEM_VECTORCALL operator*(const double2& lhs, const double rhs)
    {
        return
        {
            lhs.x * rhs,
            lhs.y * rhs
        };
    }

    GEM_INLINE double2 GEM_VECTORCALL operator/(const double2& lhs, const double rhs)
    {
        double inv = 1.0f / rhs;
        return
        {
            lhs.x * inv,
            lhs.y * inv
        };
    }

    GEM_INLINE double2 GEM_VECTORCALL operator*(const double lhs, const double2& rhs)
    {
        return
        {
            rhs.x * lhs,
            rhs.y * lhs
        };
    }

#pragma endregion

#pragma region double3

    struct double3
    {
        double x, y, z;

        double3(const double ix, const double iy, const double iz);

        double3(const double* p_components);

        double3(const double component);

        double3(const double2& component, const double z);

        double3(const double3& o);

        double3() = default;

        double length_squared() const;

        double length() const;

        double3& normalize();

        double3& safe_normalize(const double tolerance = 0.001f);

        double3& GEM_VECTORCALL operator+=(const double3& rhs);

        double3& GEM_VECTORCALL operator-=(const double3& rhs);

        double3& GEM_VECTORCALL operator*=(const double3& rhs);

        double3& GEM_VECTORCALL operator=(const double3& rhs);

        double3& operator*=(const double rhs);

        double3& operator/=(const double rhs);

        double& operator[](const unsigned int index);

        double operator[](const unsigned int index) const;

        double3 operator-() const;
    };

    double GEM_VECTORCALL dot(const double lx, const double ly, const double lz, const double rx, const double ry, const double rz);

    double GEM_VECTORCALL dot(const double3& lhs, const double3& rhs);

    double GEM_VECTORCALL length_squared(const double x, const double y, const double z);

    double GEM_VECTORCALL length_squared(const double3& val);

    double GEM_VECTORCALL length(const double x, const double y, const double z);

    double GEM_VECTORCALL length(const double3& rhs);

    double3 GEM_VECTORCALL normalize(const double3& rhs);

    double3 GEM_VECTORCALL safe_normalize(const double3& rhs, const double tolerance = 0.001f);

    double3 GEM_VECTORCALL lerp(const double3& lhs, const double3& rhs, const double t);

    double3 GEM_VECTORCALL project(const double3& lhs, const double3& rhs);

    double3 GEM_VECTORCALL rotate_axis_angle(const double3& v, const double3& a, const double angle);

    double3 GEM_VECTORCALL reject(const double3& lhs, const double3& rhs);

    double3 GEM_VECTORCALL cross(const double3& lhs, const double3& rhs);

    double3 GEM_VECTORCALL operator+(const double3& lhs, const double3& rhs);

    double3 GEM_VECTORCALL operator-(const double3& lhs, const double3& rhs);

    double3 GEM_VECTORCALL operator*(const double3& lhs, const double3& rhs);

    double3 GEM_VECTORCALL operator*(const double3& lhs, const double rhs);

    double3 GEM_VECTORCALL operator/(const double3& lhs, const double rhs);

    double3 GEM_VECTORCALL operator*(const double lhs, const double3& rhs);

    GEM_INLINE double3::double3(const double ix, const double iy, const double iz)
        : x(ix)
        , y(iy)
        , z(iz)
    {

    }

    GEM_INLINE double3::double3(const double* p_components)
        : x(p_components[0])
        , y(p_components[1])
        , z(p_components[2])
    {

    }

    GEM_INLINE double3::double3(const double component)
        : x(component)
        , y(component)
        , z(component)
    {

    }

    GEM_INLINE double3::double3(const double2& component, const double z)
        : x(component.x)
        , y(component.y)
        , z(z)
    {

    }

    GEM_INLINE double3::double3(const double3& o)
        : x(o.x)
        , y(o.y)
        , z(o.z)
    {

    }

    GEM_INLINE double double3::length_squared() const
    {
        return (x * x) + (y * y) + (z * z);
    }

    GEM_INLINE double double3::length() const
    {
        return sqrt((x * x) + (y * y) + (z * z));
    }

    GEM_INLINE double3& double3::normalize()
    {
        double il = 1.0f / sqrt((x * x) + (y * y) + (z * z));
        x *= il;
        y *= il;
        z *= il;
        return *this;
    }

    GEM_INLINE double3& double3::safe_normalize(const double tolerance /*= 0.001f*/)
    {
        double l = sqrt((x * x) + (y * y) + (z * z));
        if (l > tolerance)
        {
            double il = 1.0f / l;
            x *= il;
            y *= il;
            z *= il;
        }

        return *this;
    }

    GEM_INLINE double3& GEM_VECTORCALL double3::operator+=(const double3& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    GEM_INLINE double3& GEM_VECTORCALL double3::operator-=(const double3& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    GEM_INLINE double3& GEM_VECTORCALL double3::operator*=(const double3& rhs)
    {
        x *= rhs.x;
        y *= rhs.y;
        z *= rhs.z;
        return *this;
    }

    GEM_INLINE double3& GEM_VECTORCALL double3::operator=(const double3& rhs)
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        return *this;
    }

    GEM_INLINE double3& double3::operator*=(const double rhs)
    {
        x *= rhs;
        y *= rhs;
        z *= rhs;
        return *this;
    }

    GEM_INLINE double3& double3::operator/=(const double rhs)
    {
        double inv = 1.0f / rhs;
        x *= inv;
        y *= inv;
        z *= inv;
        return *this;
    }

    GEM_INLINE double3 double3::operator-() const
    {
        return
        {
            -x,
            -y,
            -z
        };
    }

    GEM_INLINE double& double3::operator[](const unsigned int index)
    {
        return reinterpret_cast<double*>(this)[index];
    }

    GEM_INLINE double double3::operator[](const unsigned int index) const
    {
        return reinterpret_cast<const double*>(this)[index];
    }

    GEM_INLINE double GEM_VECTORCALL dot(const double lx, const double ly, const double lz, const double rx, const double ry, const double rz)
    {
        return (lx * rx) + (ly * ry) + (lz * rz);
    }

    GEM_INLINE double GEM_VECTORCALL dot(const double3& lhs, const double3& rhs)
    {
        return (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z);
    }

    GEM_INLINE double GEM_VECTORCALL  length_squared(const double x, const double y, const double z)
    {
        return sqrt((x * x) + (y * y) + (z * z));
    }

    GEM_INLINE double GEM_VECTORCALL length_squared(const double3& val)
    {
        return (val.x * val.x) + (val.y * val.y) + (val.z * val.z);
    }

    GEM_INLINE double GEM_VECTORCALL length(const double x, const double y, const double z)
    {
        return sqrt((x * x) + (y * y) + (z * z));
    }

    GEM_INLINE double GEM_VECTORCALL length(const double3& rhs)
    {
        return sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z));
    }

    GEM_INLINE double3 GEM_VECTORCALL normalize(const double3& rhs)
    {
        double il = 1.0f / sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z));
        return
        {
            rhs.x * il,
            rhs.y * il,
            rhs.z * il
        };
    }

    GEM_INLINE double3 GEM_VECTORCALL safe_normalize(const double3& rhs, const double tolerance /*= 0.001f*/)
    {
        double l = sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z));
        double il = (l > tolerance) ? (1.0f / l) : 1.0f;
        return
        {
            rhs.x * il,
            rhs.y * il,
            rhs.z * il
        };
    }

    GEM_INLINE double3 GEM_VECTORCALL lerp(const double3& lhs, const double3& rhs, const double t)
    {
        double oneMinusT = 1.0f - t;
        return
        {
            (oneMinusT * lhs.x) + (rhs.x * t),
            (oneMinusT * lhs.y) + (rhs.y * t),
            (oneMinusT * lhs.z) + (rhs.z * t)
        };
    }

    GEM_INLINE double3 GEM_VECTORCALL project(const double3& lhs, const double3& rhs)
    {
        double d = (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z);
        double l = (rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z);
        double t = d / l;
        return
        {
            rhs.x * t,
            rhs.y * t,
            rhs.z * t
        };
    }

    GEM_INLINE double3 GEM_VECTORCALL reject(const double3& lhs, const double3& rhs)
    {
        double d = (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z);
        double l = (rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z);
        double t = d / l;
        return
        {
            lhs.x - rhs.x * t,
            lhs.y - rhs.y * t,
            lhs.z - rhs.z * t
        };
    }

    GEM_INLINE double3 GEM_VECTORCALL rotate_axis_angle(const double3& v, const double3& a, const double angle)
    {
        // v' = (v project a) + (v reject a) * cos(theta) + (a x v) * sin(theta)
        double c0 = cos(angle);
        double c1 = (1.0 - c0) * dot(v, a);
        double c2 = sin(angle);
        return (v * c0) + (a * c1) + cross(a, v) * c2;
    }

    GEM_INLINE double3 GEM_VECTORCALL cross(const double3& lhs, const double3& rhs)
    {
        return
        {
            (lhs.y * rhs.z) - (lhs.z * rhs.y),
            (lhs.z * rhs.x) - (lhs.x * rhs.z),
            (lhs.x * rhs.y) - (lhs.y * rhs.x)
        };
    }

    GEM_INLINE double3 GEM_VECTORCALL operator+(const double3& lhs, const double3& rhs)
    {
        return
        {
            lhs.x + rhs.x,
            lhs.y + rhs.y,
            lhs.z + rhs.z
        };
    }

    GEM_INLINE double3 GEM_VECTORCALL operator-(const double3& lhs, const double3& rhs)
    {
        return
        {
            lhs.x - rhs.x,
            lhs.y - rhs.y,
            lhs.z - rhs.z
        };
    }

    GEM_INLINE double3 GEM_VECTORCALL operator*(const double3& lhs, const double3& rhs)
    {
        return
        {
            lhs.x * rhs.x,
            lhs.y * rhs.y,
            lhs.z * rhs.z
        };
    }

    GEM_INLINE double3 GEM_VECTORCALL operator*(const double3& lhs, const double rhs)
    {
        return
        {
            lhs.x * rhs,
            lhs.y * rhs,
            lhs.z * rhs
        };
    }

    GEM_INLINE double3 GEM_VECTORCALL operator/(const double3& lhs, const double rhs)
    {
        double inv = 1.0f / rhs;
        return
        {
            lhs.x * inv,
            lhs.y * inv,
            lhs.z * inv
        };
    }

    GEM_INLINE double3 GEM_VECTORCALL operator*(const double lhs, const double3& rhs)
    {
        return
        {
            rhs.x * lhs,
            rhs.y * lhs,
            rhs.z * lhs
        };
    }

#pragma endregion

#pragma region double4

    struct double4
    {
        double x, y, z, w;

        double4(const double ix, const double iy, const double iz, const double iw);

        double4(const double* p_components);

        double4(const double component);

        double4(const double3& component, const double w);

        double4(const double2& component, const double z, const double w);

        double4(const double4& o);

        double4() = default;

        double length_squared() const;

        double length() const;

        double4& normalize();

        double4& safe_normalize(const double tolerance = 0.001f);

        double4& GEM_VECTORCALL operator+=(const double4& rhs);

        double4& GEM_VECTORCALL operator-=(const double4& rhs);

        double4& GEM_VECTORCALL operator*=(const double4& rhs);

        double4& GEM_VECTORCALL operator=(const double4& rhs);

        double4& operator*=(const double rhs);

        double4& operator/=(const double rhs);

        double& operator[](const unsigned int index);

        double operator[](const unsigned int index) const;

        double4 operator-() const;
    };

    double GEM_VECTORCALL dot(const double lx, const double ly, const double lz, const double lw, const double rx, const double ry, const double rz, const double rw);

    double GEM_VECTORCALL dot(const double4& lhs, const double4& rhs);

    double GEM_VECTORCALL length_squared(const double x, const double y, const double z, const double w);

    double GEM_VECTORCALL length_squared(const double4& val);

    double GEM_VECTORCALL length(const double x, const double y, const double z, const double w);

    double GEM_VECTORCALL length(const double4& rhs);

    double4 GEM_VECTORCALL normalize(const double4& rhs);

    double4 GEM_VECTORCALL safe_normalize(const double4& rhs, const double tolerance = 0.001f);

    double4 GEM_VECTORCALL lerp(const double4& lhs, const double4& rhs, const double t);

    double4 GEM_VECTORCALL project(const double4& lhs, const double4& rhs);

    double4 GEM_VECTORCALL reject(const double4& lhs, const double4& rhs);

    double4 GEM_VECTORCALL operator+(const double4& lhs, const double4& rhs);

    double4 GEM_VECTORCALL operator-(const double4& lhs, const double4& rhs);

    double4 GEM_VECTORCALL operator*(const double4& lhs, const double4& rhs);

    double4 GEM_VECTORCALL operator*(const double4& lhs, const double rhs);

    double4 GEM_VECTORCALL operator/(const double4& lhs, const double rhs);

    double4 GEM_VECTORCALL operator*(const double lhs, const double4& rhs);

    GEM_INLINE double4::double4(const double ix, const double iy, const double iz, const double iw)
        : x(ix)
        , y(iy)
        , z(iz)
        , w(iw)
    {

    }

    GEM_INLINE double4::double4(const double* p_components)
        : x(p_components[0])
        , y(p_components[1])
        , z(p_components[2])
        , w(p_components[3])
    {

    }

    GEM_INLINE double4::double4(const double component)
        : x(component)
        , y(component)
        , z(component)
        , w(component)
    {

    }

    GEM_INLINE double4::double4(const double4& o)
        : x(o.x)
        , y(o.y)
        , z(o.z)
        , w(o.w)
    {

    }

    GEM_INLINE double4::double4(const double3& component, const double w)
        : x(component.x)
        , y(component.y)
        , z(component.z)
        , w(w)
    {

    }

    GEM_INLINE double4::double4(const double2& component, const double z, const double w)
        : x(component.x)
        , y(component.y)
        , z(z)
        , w(w)
    {

    }

    GEM_INLINE double double4::length_squared() const
    {
        return (x * x) + (y * y) + (z * z) + (w * w);
    }

    GEM_INLINE double double4::length() const
    {
        return sqrt((x * x) + (y * y) + (z * z) + (w * w));
    }

    GEM_INLINE double4& double4::normalize()
    {
        double il = 1.0f / sqrt((x * x) + (y * y) + (z * z) + (w * w));
        x *= il;
        y *= il;
        z *= il;
        w *= il;
        return *this;
    }

    GEM_INLINE double4& double4::safe_normalize(const double tolerance /*= 0.001f*/)
    {
        double l = sqrt((x * x) + (y * y) + (z * z) + (w * w));
        if (l > tolerance)
        {
            double il = 1.0f / l;
            x *= il;
            y *= il;
            z *= il;
            w *= il;
        }

        return *this;
    }

    GEM_INLINE double4& GEM_VECTORCALL double4::operator+=(const double4& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        w += rhs.w;
        return *this;
    }

    GEM_INLINE double4& GEM_VECTORCALL double4::operator-=(const double4& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        w -= rhs.w;
        return *this;
    }

    GEM_INLINE double4& GEM_VECTORCALL double4::operator*=(const double4& rhs)
    {
        x *= rhs.x;
        y *= rhs.y;
        z *= rhs.z;
        w *= rhs.w;
        return *this;
    }

    GEM_INLINE double4& GEM_VECTORCALL double4::operator=(const double4& rhs)
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        w = rhs.w;
        return *this;
    }

    GEM_INLINE double4& double4::operator*=(const double rhs)
    {
        x *= rhs;
        y *= rhs;
        z *= rhs;
        w *= rhs;
        return *this;
    }

    GEM_INLINE double4& double4::operator/=(const double rhs)
    {
        double inv = 1.0f / rhs;
        x *= inv;
        y *= inv;
        z *= inv;
        w *= inv;
        return *this;
    }

    GEM_INLINE double& double4::operator[](const unsigned int index)
    {
        return reinterpret_cast<double*>(this)[index];
    }

    GEM_INLINE double double4::operator[](const unsigned int index) const
    {
        return reinterpret_cast<const double*>(this)[index];
    }

    GEM_INLINE double4 double4::operator-() const
    {
        return
        {
            -x,
            -y,
            -z,
            -w
        };
    }

    GEM_INLINE double GEM_VECTORCALL dot(const double lx, const double ly, const double lz, const double lw, const double rx, const double ry, const double rz, const double rw)
    {
        return (lx * rx) + (ly * ry) + (lz * rz) + (lw * rw);
    }

    GEM_INLINE double GEM_VECTORCALL dot(const double4& lhs, const double4& rhs)
    {
        return (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z) + (lhs.w * rhs.w);
    }

    GEM_INLINE double GEM_VECTORCALL length_squared(const double x, const double y, const double z, const double w)
    {
        return (x * x) + (y * y) + (z * z) + (w * w);
    }

    GEM_INLINE double GEM_VECTORCALL length_squared(const double4& val)
    {
        return (val.x * val.x) + (val.y * val.y) + (val.z * val.z) + (val.w * val.w);
    }

    GEM_INLINE double GEM_VECTORCALL length(const double x, const double y, const double z, const double w)
    {
        return sqrt((x * x) + (y * y) + (z * z) + (w * w));
    }

    GEM_INLINE double GEM_VECTORCALL length(const double4& rhs)
    {
        return sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w));
    }

    GEM_INLINE double4 GEM_VECTORCALL normalize(const double4& rhs)
    {
        double il = 1.0f / sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w));
        return
        {
            rhs.x * il,
            rhs.y * il,
            rhs.z * il,
            rhs.w * il
        };
    }

    GEM_INLINE double4 GEM_VECTORCALL safe_normalize(const double4& rhs, const double tolerance /*= 0.001f*/)
    {
        double l = sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w));
        double il = (l > tolerance) ? (1.0f / l) : 1.0f;
        return
        {
            rhs.x * il,
            rhs.y * il,
            rhs.z * il,
            rhs.w * il
        };
    }

    GEM_INLINE double4 GEM_VECTORCALL lerp(const double4& lhs, const double4& rhs, const double t)
    {
        double oneMinusT = 1.0f - t;
        return
        {
            (oneMinusT * lhs.x) + (rhs.x * t),
            (oneMinusT * lhs.y) + (rhs.y * t),
            (oneMinusT * lhs.z) + (rhs.z * t),
            (oneMinusT * lhs.w) + (rhs.w * t)
        };
    }

    GEM_INLINE double4 GEM_VECTORCALL project(const double4& lhs, const double4& rhs)
    {
        double d = (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z) + (lhs.w * rhs.w);
        double l = (rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w);
        double t = d / l;
        return
        {
            rhs.x * t,
            rhs.y * t,
            rhs.z * t,
            rhs.w * t
        };
    }

    GEM_INLINE double4 GEM_VECTORCALL reject(const double4& lhs, const double4& rhs)
    {
        double d = (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z) + (lhs.w * rhs.w);
        double l = (rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w);
        double t = d / l;
        return
        {
            lhs.x - rhs.x * t,
            lhs.y - rhs.y * t,
            lhs.z - rhs.z * t,
            lhs.w - rhs.w * t
        };
    }

    GEM_INLINE double4 GEM_VECTORCALL operator+(const double4& lhs, const double4& rhs)
    {
        return
        {
            lhs.x + rhs.x,
            lhs.y + rhs.y,
            lhs.z + rhs.z,
            lhs.w + rhs.w
        };
    }

    GEM_INLINE double4 GEM_VECTORCALL operator-(const double4& lhs, const double4& rhs)
    {
        return
        {
            lhs.x - rhs.x,
            lhs.y - rhs.y,
            lhs.z - rhs.z,
            lhs.w - rhs.w
        };
    }

    GEM_INLINE double4 GEM_VECTORCALL operator*(const double4& lhs, const double4& rhs)
    {
        return
        {
            lhs.x * rhs.x,
            lhs.y * rhs.y,
            lhs.z * rhs.z,
            lhs.w * rhs.w
        };
    }

    GEM_INLINE double4 GEM_VECTORCALL operator*(const double4& lhs, const double rhs)
    {
        return
        {
            lhs.x * rhs,
            lhs.y * rhs,
            lhs.z * rhs,
            lhs.w * rhs
        };
    }

    GEM_INLINE double4 GEM_VECTORCALL operator/(const double4& lhs, const double rhs)
    {
        double inv = 1.0f / rhs;
        return
        {
            lhs.x * inv,
            lhs.y * inv,
            lhs.z * inv,
            lhs.w * inv
        };
    }

    GEM_INLINE double4 GEM_VECTORCALL operator*(const double lhs, const double4& rhs)
    {
        return
        {
            rhs.x * lhs,
            rhs.y * lhs,
            rhs.z * lhs,
            rhs.w * lhs
        };
    }

#pragma endregion 

#pragma region int2

    struct int2
    {
        int x, y;

        int2(const int ix, const int iy);

        int2(const int* p_components);

        int2(const int component);

        int2(const int2& o);

        int2() = default;

        int length() const;

        int2& GEM_VECTORCALL operator+=(const int2& rhs);

        int2& GEM_VECTORCALL operator-=(const int2& rhs);

        int2& GEM_VECTORCALL operator*=(const int2& rhs);

        int2& GEM_VECTORCALL operator=(const int2& rhs);

        int2& operator*=(const int rhs);

        int2& operator/=(const int rhs);

        int& operator[](const unsigned int index);

        int operator[](const unsigned int index) const;

        int2 operator-() const;
    };

    int GEM_VECTORCALL dot(const int2& lhs, const int2& rhs);

    int GEM_VECTORCALL length(const int2& rhs);

    int2 GEM_VECTORCALL operator+(const int2& lhs, const int2& rhs);

    int2 GEM_VECTORCALL operator-(const int2& lhs, const int2& rhs);

    int2 GEM_VECTORCALL operator*(const int2& lhs, const int2& rhs);

    int2 GEM_VECTORCALL operator*(const int2& lhs, const int rhs);

    int2 GEM_VECTORCALL operator/(const int2& lhs, const int rhs);

    int2 GEM_VECTORCALL operator*(const int lhs, const int2& rhs);

    GEM_INLINE int2::int2(const int ix, const int iy)
        : x(ix)
        , y(iy)
    {

    }

    GEM_INLINE int2::int2(const int* p_components)
        : x(p_components[0])
        , y(p_components[1])
    {

    }

    GEM_INLINE int2::int2(const int component)
        : x(component)
        , y(component)
    {

    }

    GEM_INLINE int2::int2(const int2& o)
        : x(o.x)
        , y(o.y)
    {

    }

    GEM_INLINE int int2::length() const
    {
        return abs(x) + abs(y);
    }

    GEM_INLINE int2& GEM_VECTORCALL int2::operator+=(const int2& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    GEM_INLINE int2& GEM_VECTORCALL int2::operator-=(const int2& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    GEM_INLINE int2& GEM_VECTORCALL int2::operator*=(const int2& rhs)
    {
        x *= rhs.x;
        y *= rhs.y;
        return *this;
    }

    GEM_INLINE int2& GEM_VECTORCALL int2::operator=(const int2& rhs)
    {
        x = rhs.x;
        y = rhs.y;
        return *this;
    }

    GEM_INLINE int2& int2::operator*=(const int rhs)
    {
        x *= rhs;
        y *= rhs;
        return *this;
    }

    GEM_INLINE int2& int2::operator/=(const int rhs)
    {
        int inv = 1.0f / rhs;
        x *= inv;
        y *= inv;
        return *this;
    }

    GEM_INLINE int& int2::operator[](const unsigned int index)
    {
        return reinterpret_cast<int*>(this)[index];
    }

    GEM_INLINE int int2::operator[](const unsigned int index) const
    {
        return reinterpret_cast<const int*>(this)[index];
    }

    GEM_INLINE int2 int2::operator-() const
    {
        return
        {
            -x,
            -y
        };
    }

    GEM_INLINE int GEM_VECTORCALL dot(const int2& lhs, const int2& rhs)
    {
        return (lhs.x * rhs.x) + (lhs.y * rhs.y);
    }

    GEM_INLINE int GEM_VECTORCALL length(const int2& rhs)
    {
        return abs(rhs.x) + abs(rhs.y);
    }

    GEM_INLINE int2 GEM_VECTORCALL lerp(const int2& lhs, const int2& rhs, const int t)
    {
        int oneMinusT = 1.0f - t;
        return
        {
            (oneMinusT * lhs.x) + (rhs.x * t),
            (oneMinusT * lhs.y) + (rhs.y * t)
        };
    }

    GEM_INLINE int2 GEM_VECTORCALL operator+(const int2& lhs, const int2& rhs)
    {
        return
        {
            lhs.x + rhs.x,
            lhs.y + rhs.y
        };
    }

    GEM_INLINE int2 GEM_VECTORCALL operator-(const int2& lhs, const int2& rhs)
    {
        return
        {
            lhs.x - rhs.x,
            lhs.y - rhs.y
        };
    }

    GEM_INLINE int2 GEM_VECTORCALL operator*(const int2& lhs, const int2& rhs)
    {
        return
        {
            lhs.x * rhs.x,
            lhs.y * rhs.y
        };
    }

    GEM_INLINE int2 GEM_VECTORCALL operator*(const int2& lhs, const int rhs)
    {
        return
        {
            lhs.x * rhs,
            lhs.y * rhs
        };
    }

    GEM_INLINE int2 GEM_VECTORCALL operator/(const int2& lhs, const int rhs)
    {
        int inv = 1.0f / rhs;
        return
        {
            lhs.x * inv,
            lhs.y * inv
        };
    }

    GEM_INLINE int2 GEM_VECTORCALL operator*(const int lhs, const int2& rhs)
    {
        return
        {
            rhs.x * lhs,
            rhs.y * lhs
        };
    }

#pragma endregion

#pragma region int3

    struct int3
    {
        int x, y, z;

        int3(const int ix, const int iy, const int iz);

        int3(const int* p_components);

        int3(const int component);

        int3(const int2& component, const int z);

        int3(const int3& o);

        int3() = default;

        int length() const;

        int3& GEM_VECTORCALL operator+=(const int3& rhs);

        int3& GEM_VECTORCALL operator-=(const int3& rhs);

        int3& GEM_VECTORCALL operator*=(const int3& rhs);

        int3& GEM_VECTORCALL operator=(const int3& rhs);

        int3& operator*=(const int rhs);

        int3& operator/=(const int rhs);

        int& operator[](const unsigned int index);

        int operator[](const unsigned int index) const;

        int3 operator-() const;
    };

    int GEM_VECTORCALL dot(const int3& lhs, const int3& rhs);

    int GEM_VECTORCALL length(const int3& rhs);

    int3 GEM_VECTORCALL cross(const int3& lhs, const int3& rhs);

    int3 GEM_VECTORCALL operator+(const int3& lhs, const int3& rhs);

    int3 GEM_VECTORCALL operator-(const int3& lhs, const int3& rhs);

    int3 GEM_VECTORCALL operator*(const int3& lhs, const int3& rhs);

    int3 GEM_VECTORCALL operator*(const int3& lhs, const int rhs);

    int3 GEM_VECTORCALL operator/(const int3& lhs, const int rhs);

    int3 GEM_VECTORCALL operator*(const int lhs, const int3& rhs);

    GEM_INLINE int3::int3(const int ix, const int iy, const int iz)
        : x(ix)
        , y(iy)
        , z(iz)
    {

    }

    GEM_INLINE int3::int3(const int* p_components)
        : x(p_components[0])
        , y(p_components[1])
        , z(p_components[2])
    {

    }

    GEM_INLINE int3::int3(const int component)
        : x(component)
        , y(component)
        , z(component)
    {

    }

    GEM_INLINE int3::int3(const int2& component, const int z)
        : x(component.x)
        , y(component.y)
        , z(z)
    {

    }

    GEM_INLINE int3::int3(const int3& o)
        : x(o.x)
        , y(o.y)
        , z(o.z)
    {

    }

    GEM_INLINE int int3::length() const
    {
        return abs(x) + abs(y) + abs(z);
    }

    GEM_INLINE int3& GEM_VECTORCALL int3::operator+=(const int3& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    GEM_INLINE int3& GEM_VECTORCALL int3::operator-=(const int3& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    GEM_INLINE int3& GEM_VECTORCALL int3::operator*=(const int3& rhs)
    {
        x *= rhs.x;
        y *= rhs.y;
        z *= rhs.z;
        return *this;
    }

    GEM_INLINE int3& GEM_VECTORCALL int3::operator=(const int3& rhs)
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        return *this;
    }

    GEM_INLINE int3& int3::operator*=(const int rhs)
    {
        x *= rhs;
        y *= rhs;
        z *= rhs;
        return *this;
    }

    GEM_INLINE int3& int3::operator/=(const int rhs)
    {
        int inv = 1.0f / rhs;
        x *= inv;
        y *= inv;
        z *= inv;
        return *this;
    }

    GEM_INLINE int3 int3::operator-() const
    {
        return
        {
            -x,
            -y,
            -z
        };
    }

    GEM_INLINE int& int3::operator[](const unsigned int index)
    {
        return reinterpret_cast<int*>(this)[index];
    }

    GEM_INLINE int int3::operator[](const unsigned int index) const
    {
        return reinterpret_cast<const int*>(this)[index];
    }

    GEM_INLINE int GEM_VECTORCALL dot(const int3& lhs, const int3& rhs)
    {
        return (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z);
    }

    GEM_INLINE int GEM_VECTORCALL length(const int3& rhs)
    {
        return abs(rhs.x) + abs(rhs.y) + abs(rhs.z);
    }

    GEM_INLINE int3 GEM_VECTORCALL cross(const int3& lhs, const int3& rhs)
    {
        return
        {
            (lhs.y * rhs.z) - (lhs.z * rhs.y),
            (lhs.z * rhs.x) - (lhs.x * rhs.z),
            (lhs.x * rhs.y) - (lhs.y * rhs.x)
        };
    }

    GEM_INLINE int3 GEM_VECTORCALL operator+(const int3& lhs, const int3& rhs)
    {
        return
        {
            lhs.x + rhs.x,
            lhs.y + rhs.y,
            lhs.z + rhs.z
        };
    }

    GEM_INLINE int3 GEM_VECTORCALL operator-(const int3& lhs, const int3& rhs)
    {
        return
        {
            lhs.x - rhs.x,
            lhs.y - rhs.y,
            lhs.z - rhs.z
        };
    }

    GEM_INLINE int3 GEM_VECTORCALL operator*(const int3& lhs, const int3& rhs)
    {
        return
        {
            lhs.x * rhs.x,
            lhs.y * rhs.y,
            lhs.z * rhs.z
        };
    }

    GEM_INLINE int3 GEM_VECTORCALL operator*(const int3& lhs, const int rhs)
    {
        return
        {
            lhs.x * rhs,
            lhs.y * rhs,
            lhs.z * rhs
        };
    }

    GEM_INLINE int3 GEM_VECTORCALL operator/(const int3& lhs, const int rhs)
    {
        int inv = 1.0f / rhs;
        return
        {
            lhs.x * inv,
            lhs.y * inv,
            lhs.z * inv
        };
    }

    GEM_INLINE int3 GEM_VECTORCALL operator*(const int lhs, const int3& rhs)
    {
        return
        {
            rhs.x * lhs,
            rhs.y * lhs,
            rhs.z * lhs
        };
    }

#pragma endregion

#pragma region int4

    struct int4
    {
        int x, y, z, w;

        int4(const int ix, const int iy, const int iz, const int iw);

        int4(const int* p_components);

        int4(const int component);

        int4(const int3& component, const int w);

        int4(const int2& component, const int z, const int w);

        int4(const int4& o);

        int4() = default;

        int length() const;

        int4& GEM_VECTORCALL operator+=(const int4& rhs);

        int4& GEM_VECTORCALL operator-=(const int4& rhs);

        int4& GEM_VECTORCALL operator*=(const int4& rhs);

        int4& GEM_VECTORCALL operator=(const int4& rhs);

        int4& operator*=(const int rhs);

        int4& operator/=(const int rhs);

        int& operator[](const unsigned int index);

        int operator[](const unsigned int index) const;

        int4 operator-() const;
    };

    int GEM_VECTORCALL dot(const int4& lhs, const int4& rhs);

    int GEM_VECTORCALL length(const int4& rhs);

    int4 GEM_VECTORCALL operator+(const int4& lhs, const int4& rhs);

    int4 GEM_VECTORCALL operator-(const int4& lhs, const int4& rhs);

    int4 GEM_VECTORCALL operator*(const int4& lhs, const int4& rhs);

    int4 GEM_VECTORCALL operator*(const int4& lhs, const int rhs);

    int4 GEM_VECTORCALL operator/(const int4& lhs, const int rhs);

    int4 GEM_VECTORCALL operator*(const int lhs, const int4& rhs);

    GEM_INLINE int4::int4(const int ix, const int iy, const int iz, const int iw)
        : x(ix)
        , y(iy)
        , z(iz)
        , w(iw)
    {

    }

    GEM_INLINE int4::int4(const int* p_components)
        : x(p_components[0])
        , y(p_components[1])
        , z(p_components[2])
        , w(p_components[3])
    {

    }

    GEM_INLINE int4::int4(const int component)
        : x(component)
        , y(component)
        , z(component)
        , w(component)
    {

    }

    GEM_INLINE int4::int4(const int4& o)
        : x(o.x)
        , y(o.y)
        , z(o.z)
        , w(o.w)
    {

    }

    GEM_INLINE int4::int4(const int3& component, const int w)
        : x(component.x)
        , y(component.y)
        , z(component.z)
        , w(w)
    {

    }

    GEM_INLINE int4::int4(const int2& component, const int z, const int w)
        : x(component.x)
        , y(component.y)
        , z(z)
        , w(w)
    {

    }

    GEM_INLINE int int4::length() const
    {
        return abs(x) + abs(y) + abs(z) + abs(w);
    }

    GEM_INLINE int4& GEM_VECTORCALL int4::operator+=(const int4& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        w += rhs.w;
        return *this;
    }

    GEM_INLINE int4& GEM_VECTORCALL int4::operator-=(const int4& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        w -= rhs.w;
        return *this;
    }

    GEM_INLINE int4& GEM_VECTORCALL int4::operator*=(const int4& rhs)
    {
        x *= rhs.x;
        y *= rhs.y;
        z *= rhs.z;
        w *= rhs.w;
        return *this;
    }

    GEM_INLINE int4& GEM_VECTORCALL int4::operator=(const int4& rhs)
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        w = rhs.w;
        return *this;
    }

    GEM_INLINE int4& int4::operator*=(const int rhs)
    {
        x *= rhs;
        y *= rhs;
        z *= rhs;
        w *= rhs;
        return *this;
    }

    GEM_INLINE int4& int4::operator/=(const int rhs)
    {
        int inv = 1.0f / rhs;
        x *= inv;
        y *= inv;
        z *= inv;
        w *= inv;
        return *this;
    }

    GEM_INLINE int& int4::operator[](const unsigned int index)
    {
        return reinterpret_cast<int*>(this)[index];
    }

    GEM_INLINE int int4::operator[](const unsigned int index) const
    {
        return reinterpret_cast<const int*>(this)[index];
    }

    GEM_INLINE int4 int4::operator-() const
    {
        return
        {
            -x,
            -y,
            -z,
            -w
        };
    }

    GEM_INLINE int GEM_VECTORCALL dot(const int4& lhs, const int4& rhs)
    {
        return (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z) + (lhs.w * rhs.w);
    }

    GEM_INLINE int GEM_VECTORCALL length(const int4& rhs)
    {
        return abs(rhs.x) + abs(rhs.y) + abs(rhs.z) + abs(rhs.w);
    }

    GEM_INLINE int4 GEM_VECTORCALL operator+(const int4& lhs, const int4& rhs)
    {
        return
        {
            lhs.x + rhs.x,
            lhs.y + rhs.y,
            lhs.z + rhs.z,
            lhs.w + rhs.w
        };
    }

    GEM_INLINE int4 GEM_VECTORCALL operator-(const int4& lhs, const int4& rhs)
    {
        return
        {
            lhs.x - rhs.x,
            lhs.y - rhs.y,
            lhs.z - rhs.z,
            lhs.w - rhs.w
        };
    }

    GEM_INLINE int4 GEM_VECTORCALL operator*(const int4& lhs, const int4& rhs)
    {
        return
        {
            lhs.x * rhs.x,
            lhs.y * rhs.y,
            lhs.z * rhs.z,
            lhs.w * rhs.w
        };
    }

    GEM_INLINE int4 GEM_VECTORCALL operator*(const int4& lhs, const int rhs)
    {
        return
        {
            lhs.x * rhs,
            lhs.y * rhs,
            lhs.z * rhs,
            lhs.w * rhs
        };
    }

    GEM_INLINE int4 GEM_VECTORCALL operator/(const int4& lhs, const int rhs)
    {
        int inv = 1.0f / rhs;
        return
        {
            lhs.x * inv,
            lhs.y * inv,
            lhs.z * inv,
            lhs.w * inv
        };
    }

    GEM_INLINE int4 GEM_VECTORCALL operator*(const int lhs, const int4& rhs)
    {
        return
        {
            rhs.x * lhs,
            rhs.y * lhs,
            rhs.z * lhs,
            rhs.w * lhs
        };
    }

#pragma endregion 
}
