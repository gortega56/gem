#pragma once
#include "matrix.h"

namespace gem
{
#pragma region quatf

    struct quatf
    {
        float x, y, z, w;

        static quatf identity();

        static quatf GEM_VECTORCALL rotate_euler(const float pitch, const float yaw, const float roll);

        static quatf GEM_VECTORCALL rotate_axis_angle(const float3& axis, const float angle);

        static quatf GEM_VECTORCALL rotate_matrix(const float3x3& m);

        static quatf GEM_VECTORCALL rotate_align(const float3 u, const float3 v);

        quatf(const float ix, const float iy, const float iz, const float iw);

        quatf(const float* o);

        quatf(const quatf& o);

        quatf() = default;

        float length_squared() const;

        float length() const;

        quatf& normalize();

        quatf& safe_normalize(const float tolerance = 0.001f);

        quatf& conjugate();

        quatf& inverse();

        float3x3 matrix3x3() const;

        float4x3 matrix4x3() const;

        float4x4 matrix4x4() const;

        float3 euler() const;

        float4 axis_angle() const;

        float3 axis_x() const;

        float3 axis_y() const;

        float3 axis_z() const;

        float3 GEM_VECTORCALL transform_point(const float3& p) const;

        quatf& GEM_VECTORCALL operator=(const quatf& rhs);

        quatf& GEM_VECTORCALL operator+=(const quatf& rhs);

        quatf& GEM_VECTORCALL operator-=(const quatf& rhs);

        quatf& GEM_VECTORCALL operator*=(const quatf& rhs);

        quatf& GEM_VECTORCALL operator*=(const float rhs);

        quatf& GEM_VECTORCALL operator/=(const float rhs);

        quatf operator-() const;
    };

    float GEM_VECTORCALL dot(const quatf& lhs, const quatf& rhs);

    float GEM_VECTORCALL length(const quatf& rhs);

    float GEM_VECTORCALL length_squared(const quatf& val);

    quatf GEM_VECTORCALL normalize(const quatf& rhs);

    quatf GEM_VECTORCALL safe_normalize(const quatf& rhs, const float tolerance = 0.001f);

    quatf GEM_VECTORCALL conjugate(const quatf& rhs);

    quatf GEM_VECTORCALL inverse(const quatf& rhs);

    quatf GEM_VECTORCALL slerp(quatf q0, quatf q1, const float t);

    quatf GEM_VECTORCALL operator+(const quatf& lhs, const quatf& rhs);

    quatf GEM_VECTORCALL operator-(const quatf& lhs, const quatf& rhs);

    quatf GEM_VECTORCALL operator*(const quatf& lhs, const quatf& rhs);

    quatf GEM_VECTORCALL operator*(const quatf& lhs, const float rhs);

    quatf GEM_VECTORCALL operator*(const float lhs, const quatf& rhs);

    quatf GEM_VECTORCALL operator/(const quatf& lhs, const float rhs);

    GEM_INLINE quatf quatf::identity()
    {
        return { 0, 0, 0, 1 };
    }

    GEM_INLINE quatf GEM_VECTORCALL quatf::rotate_euler(const float pitch, const float yaw, const float roll)
    {
        float halfRoll = roll * 0.5f;
        float halfPitch = pitch * 0.5f;
        float halfYaw = yaw * 0.5f;
        float cosHalfRoll = cosf(halfRoll);
        float cosHalfPitch = cosf(halfPitch);
        float cosHalfYaw = cosf(halfYaw);
        float sinHalfRoll = sinf(halfRoll);
        float sinHalfPitch = sinf(halfPitch);
        float sinHalfYaw = sinf(halfYaw);

        return
        {
            (cosHalfYaw * sinHalfPitch * cosHalfRoll) + (sinHalfYaw * cosHalfPitch * sinHalfRoll),
            (sinHalfYaw * cosHalfPitch * cosHalfRoll) - (cosHalfYaw * sinHalfPitch * sinHalfRoll),
            (cosHalfYaw * cosHalfPitch * sinHalfRoll) - (sinHalfYaw * sinHalfPitch * cosHalfRoll),
            (cosHalfYaw * cosHalfPitch * cosHalfRoll) + (sinHalfYaw * sinHalfPitch * sinHalfRoll)
        };
    }

    GEM_INLINE quatf GEM_VECTORCALL quatf::rotate_axis_angle(const float3& axis, const float angle)
    {
        float a = angle * 0.5f;
        float s = sinf(a);
        float c = cosf(a);

        return
        {
            axis.x * s,
            axis.y * s,
            axis.z * s,
            c
        };
    }

    GEM_INLINE quatf GEM_VECTORCALL quatf::rotate_matrix(const float3x3& m)
    {
        quatf o;
        const float m11 = m.m00, m12 = m.m01, m13 = m.m02;
        const float m21 = m.m10, m22 = m.m11, m23 = m.m12;
        const float m31 = m.m20, m32 = m.m21, m33 = m.m22;

        // Determine which of w, x, y or z has the largest absolute value
        float fourWSquaredMinus1 = +m11 + m22 + m33;
        float fourXSquaredMinus1 = +m11 - m22 - m33;
        float fourYSquaredMinus1 = -m11 + m22 - m33;
        float fourZSquaredMinus1 = -m11 - m22 + m33;

        int biggestIndex = 0;
        float fourBiggestSquardeMinus1 = fourWSquaredMinus1;
        if (fourXSquaredMinus1 > fourBiggestSquardeMinus1)
        {
            fourBiggestSquardeMinus1 = fourXSquaredMinus1;
            biggestIndex = 1;
        }
        if (fourYSquaredMinus1 > fourBiggestSquardeMinus1)
        {
            fourBiggestSquardeMinus1 = fourYSquaredMinus1;
            biggestIndex = 2;
        }
        if (fourZSquaredMinus1 > fourBiggestSquardeMinus1)
        {
            fourBiggestSquardeMinus1 = fourZSquaredMinus1;
            biggestIndex = 3;
        }

        float biggestVal = sqrt(fourBiggestSquardeMinus1 + 1) * .5f;
        float mult = 0.25f / biggestVal;

        switch (biggestIndex)
        {
        case 0:
        {
            o.x = (m23 - m32) * mult;
            o.y = (m31 - m13) * mult;
            o.z = (m12 - m21) * mult;
            o.w = biggestVal;
            break;
        }
        case 1:
        {
            o.x = biggestVal;
            o.y = (m12 + m21) * mult;
            o.z = (m31 + m13) * mult;
            o.w = (m23 - m32) * mult;
            break;
        }
        case 2:
        {
            o.x = (m12 + m21) * mult;
            o.y = biggestVal;
            o.z = (m23 + m32) * mult;
            o.w = (m31 - m13) * mult;
            break;
        }
        case 3:
        {
            o.x = (m31 + m13) * mult;
            o.y = (m23 + m32) * mult;
            o.z = biggestVal;
            o.w = (m12 - m21) * mult;
            break;
        }
        default:
            o.x = 0.0f;
            o.y = 0.0f;
            o.z = 0.0f;
            o.w = 1.0f;
            break;
        }
        return o;
    }

    GEM_INLINE quatf GEM_VECTORCALL quatf::rotate_align(const float3 u, const float3 v)
    {
        float c = sqrtf((1.f + dot(u, v)) * 0.5f);
        float s = 1.f / (2.f * c);
        float3 a = cross(u, v);
        return
        {
            a.x * s,
            a.y * s,
            a.z * s,
            c
        };
    }

    GEM_INLINE quatf::quatf(const float ix, const float iy, const float iz, const float iw)
        : x(ix)
        , y(iy)
        , z(iz)
        , w(iw)
    {

    }

    GEM_INLINE quatf::quatf(const float* o)
        : x(o[0])
        , y(o[1])
        , z(o[2])
        , w(o[3])
    {

    }

    GEM_INLINE quatf::quatf(const quatf& o)
        : x(o.x)
        , y(o.y)
        , z(o.z)
        , w(o.w)
    {

    }

    GEM_INLINE float quatf::length_squared() const
    {
        return (x * x) + (y * y) + (z * z) + (w * w);
    }

    GEM_INLINE float quatf::length() const
    {
        return sqrtf((x * x) + (y * y) + (z * z) + (w * w));
    }

    GEM_INLINE quatf& quatf::normalize()
    {
        float il = 1.0f / sqrtf((x * x) + (y * y) + (z * z) + (w * w));
        x *= il;
        y *= il;
        z *= il;
        w *= il;
        return *this;
    }

    GEM_INLINE quatf& quatf::safe_normalize(const float tolerance /*= 0.001f*/)
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

    GEM_INLINE quatf& quatf::conjugate()
    {
        x = -x;
        y = -y;
        z = -z;
        return *this;
    }

    GEM_INLINE quatf& quatf::inverse()
    {
        float il = 1.0f / sqrtf((x * x) + (y * y) + (z * z) + (w * w));
        x = -x * il;
        y = -y * il;
        z = -z * il;
        w = w * il;
        return *this;
    }

    GEM_INLINE float3x3 quatf::matrix3x3() const
    {
        float x2 = x * x;
        float y2 = y * y;
        float z2 = z * z;
        float wx = w * x;
        float wy = w * y;
        float wz = w * z;
        float xy = x * y;
        float xz = x * z;
        float yz = y * z;
        return
        {
            1.0f - (2.0f * y2) - (2.0f * z2),
            (2.0f * xy) + (2.0f * wz),
            (2.0f * xz) - (2.0f * wy),
            (2.0f * xy) - (2.0f * wz),
            1.0f - (2.0f * x2) - (2.0f * z2),
            (2.0f * yz) + (2.0f * wx),
            (2.0f * xz) + (2.0f * wy),
            (2.0f * yz) - (2.0f * wx),
            1.0f - (2.0f * x2) - (2.0f * y2)
        };
    }

    GEM_INLINE float4x3 quatf::matrix4x3() const
    {
        float x2 = x * x;
        float y2 = y * y;
        float z2 = z * z;
        float wx = w * x;
        float wy = w * y;
        float wz = w * z;
        float xy = x * y;
        float xz = x * z;
        float yz = y * z;

        return
        {
            1.0f - (2.0f * y2) - (2.0f * z2),
            (2.0f * xy) + (2.0f * wz),
            (2.0f * xz) - (2.0f * wy),
            (2.0f * xy) - (2.0f * wz),
            1.0f - (2.0f * x2) - (2.0f * z2),
            (2.0f * yz) + (2.0f * wx),
            (2.0f * xz) + (2.0f * wy),
            (2.0f * yz) - (2.0f * wx),
            1.0f - (2.0f * x2) - (2.0f * y2),
            0.0f,
            0.0f,
            0.0f
        };
    }

    GEM_INLINE float4x4 quatf::matrix4x4() const
    {
        float x2 = x * x;
        float y2 = y * y;
        float z2 = z * z;
        float wx = w * x;
        float wy = w * y;
        float wz = w * z;
        float xy = x * y;
        float xz = x * z;
        float yz = y * z;

        return
        {
            1.0f - (2.0f * y2) - (2.0f * z2),
            (2.0f * xy) + (2.0f * wz),
            (2.0f * xz) - (2.0f * wy),
            0.0f,
            (2.0f * xy) - (2.0f * wz),
            1.0f - (2.0f * x2) - (2.0f * z2),
            (2.0f * yz) + (2.0f * wx),
            0.0f,
            (2.0f * xz) + (2.0f * wy),
            (2.0f * yz) - (2.0f * wx),
            1.0f - (2.0f * x2) - (2.0f * y2),
            0.0f,
            0.0f,
            0.0f,
            0.0f,
            1.0f
        };
    }

    GEM_INLINE float3 quatf::euler() const
    {
        float3 o;
        auto sp = -2.0f * ((y * z) - (w * x));
        if (fabs(sp) > 0.9999f)
        {
            o.x = 1.570796f * sp;
            o.y = atan2((-x * z) + (w * y), 0.5f - (y * y) - (z * z));
            o.z = 0.0f;
        }
        else
        {
            o.x = asin(sp);
            o.y = atan2((x * z) + (w * y), 0.5f - (x * x) - (y * y));
            o.z = atan2((x * y) + (w * z), 0.5f - (x * x) - (z * z));
        }
        return o;
    }

    GEM_INLINE float4 quatf::axis_angle() const
    {
        float half_angle = acosf(w);
        float s = 1.f / sinf(half_angle);
        return
        {
            x * s,
            y * s,
            z * s,
            w * half_angle * 2.f
        };
    }

    GEM_INLINE float3 quatf::axis_x() const
    {
        float y2 = y * y;
        float z2 = z * z;
        float wy = w * y;
        float wz = w * z;
        float xy = x * y;
        float xz = x * z;
        return  
        {
            1.0f - (2.0f * y2) - (2.0f * z2),
            (2.0f * xy) + (2.0f * wz),
            (2.0f * xz) - (2.0f * wy)
        };
    }

    GEM_INLINE float3 quatf::axis_y() const
    {
        float x2 = x * x;
        float z2 = z * z;
        float wx = w * x;
        float wz = w * z;
        float xy = x * y;
        float yz = y * z;
        return
        {
            (2.0f * xy) - (2.0f * wz),
            1.0f - (2.0f * x2) - (2.0f * z2),
            (2.0f * yz) + (2.0f * wx),
        };
    }

    GEM_INLINE float3 quatf::axis_z() const
    {
        float x2 = x * x;
        float y2 = y * y;
        float wx = w * x;
        float wy = w * y;
        float xz = x * z;
        float yz = y * z;
        return
        {
            (2.0f * xz) + (2.0f * wy),
            (2.0f * yz) - (2.0f * wx),
            1.0f - (2.0f * x2) - (2.0f * y2)
        };
    }

    GEM_INLINE float3 GEM_VECTORCALL quatf::transform_point(const float3& p) const
    {
        // Given q is a unit vector:
        // qvq* = (b + c)(v + 0)(-b + c)
        //      = v + 2c(bxv) + 2(bxbxv)
        float3 v = { x, y, z };
        float3 bxv = cross(v, p);
        float3 bxbxv = cross(v, bxv);
        return 
        {
            p.x + ((bxv.x * w) + bxbxv.x) * 2.0f,
            p.y + ((bxv.y * w) + bxbxv.y) * 2.0f,
            p.z + ((bxv.z * w) + bxbxv.z) * 2.0f
        };
    }

    GEM_INLINE quatf& GEM_VECTORCALL quatf::operator=(const quatf& rhs)
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        w = rhs.w;
        return *this;
    }

    GEM_INLINE quatf& GEM_VECTORCALL quatf::operator+=(const quatf& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        w += rhs.w;
        return *this;
    }

    GEM_INLINE quatf& GEM_VECTORCALL quatf::operator-=(const quatf& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        w -= rhs.w;
        return *this;
    }

    GEM_INLINE quatf& GEM_VECTORCALL quatf::operator*=(const quatf& rhs)
    {
        float ow = (w * rhs.w) - ((x * rhs.x) + (y * rhs.y) + (z * rhs.z));
        float ox = (w * rhs.x) + (rhs.w * x) + ((y * rhs.z) - (z * rhs.y));
        float oy = (w * rhs.y) + (rhs.w * y) + ((z * rhs.x) - (x * rhs.z));
        float oz = (w * rhs.z) + (rhs.w * z) + ((x * rhs.y) - (y * rhs.x));
        x = ox;
        y = oy;
        z = oz;
        w = ow;
        return *this;
    }

    GEM_INLINE quatf& GEM_VECTORCALL quatf::operator*=(const float rhs)
    {
        x *= rhs;
        y *= rhs;
        z *= rhs;
        w *= rhs;
        return *this;
    }

    GEM_INLINE quatf& GEM_VECTORCALL quatf::operator/=(const float rhs)
    {
        float inv = 1.0f / rhs;
        x *= inv;
        y *= inv;
        z *= inv;
        w *= inv;
        return *this;
    }

    GEM_INLINE quatf quatf::operator-() const
    {
        return
        {
            -x,
            -y,
            -z,
            -w
        };
    }

    GEM_INLINE float GEM_VECTORCALL dot(const quatf& lhs, const quatf& rhs)
    {
        return (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z) + (lhs.w * rhs.w);
    }

    GEM_INLINE float GEM_VECTORCALL length_squared(const quatf& val)
    {
        return (val.x * val.x) + (val.y * val.y) + (val.z * val.z) + (val.w * val.w);
    }

    GEM_INLINE float GEM_VECTORCALL length(const quatf& rhs)
    {
        return sqrtf((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w));
    }

    GEM_INLINE quatf GEM_VECTORCALL normalize(const quatf& rhs)
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

    GEM_INLINE quatf GEM_VECTORCALL safe_normalize(const quatf& rhs, const float tolerance /*= 0.001f*/)
    {
        float l = sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w));
        float il = (l > tolerance) ? (1.0f / l) : 1.0f;
        return
        {
            rhs.x * il,
            rhs.y * il,
            rhs.z * il,
            rhs.w * il
        };
    }

    GEM_INLINE quatf GEM_VECTORCALL conjugate(const quatf& rhs)
    {
        return { -rhs.x, -rhs.y, -rhs.z, rhs.w };
    }

    GEM_INLINE quatf GEM_VECTORCALL inverse(const quatf& rhs)
    {
        float il = 1.0f / sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w));
        return { -rhs.x * il, -rhs.y * il, -rhs.z * il, rhs.w * il };
    }

    GEM_INLINE quatf GEM_VECTORCALL slerp(quatf q0, quatf q1, const float t)
    {
        float cosAngle = (q0.x * q1.x) + (q0.y * q1.y) + (q0.z * q1.z) + (q0.w * q1.w);
        if (cosAngle < 0.0f) {
            q1 = -q1;
            cosAngle = -cosAngle;
        }

        float k0, k1;

        // Check for divide by zero
        if (cosAngle > 0.9999f) {
            k0 = 1.0f - t;
            k1 = t;
        }
        else {
            float angle = acosf(cosAngle);
            float oneOverSinAngle = 1.0f / sinf(angle);

            k0 = ((sinf(1.0f - t) * angle) * oneOverSinAngle);
            k1 = (sinf(t * angle) * oneOverSinAngle);
        }

        q0 = q0 * k0;
        q1 = q1 * k1;

        return q0 + q1;
    }

    GEM_INLINE quatf GEM_VECTORCALL operator+(const quatf& lhs, const quatf& rhs)
    {
        return
        {
            lhs.x + rhs.x,
            lhs.y + rhs.y,
            lhs.z + rhs.z,
            lhs.w + rhs.w
        };
    }

    GEM_INLINE quatf GEM_VECTORCALL operator-(const quatf& lhs, const quatf& rhs)
    {
        return
        {
            lhs.x - rhs.x,
            lhs.y - rhs.y,
            lhs.z - rhs.z,
            lhs.w - rhs.w
        };
    }

    GEM_INLINE quatf GEM_VECTORCALL operator*(const quatf& lhs, const quatf& rhs)
    {
        return
        {
            (lhs.w * rhs.x) + (lhs.x * rhs.w) + ((lhs.y * rhs.z) - (lhs.z * rhs.y)),
            (lhs.w * rhs.y) + (lhs.y * rhs.w) + ((lhs.z * rhs.x) - (lhs.x * rhs.z)),
            (lhs.w * rhs.z) + (lhs.z * rhs.w) + ((lhs.x * rhs.y) - (lhs.y * rhs.x)),
            (lhs.w * rhs.w) - ((lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z))
        };
    }

    GEM_INLINE quatf GEM_VECTORCALL operator*(const quatf& lhs, const float rhs)
    {
        return
        {
            lhs.x * rhs,
            lhs.y * rhs,
            lhs.z * rhs,
            lhs.w * rhs
        };
    }

    GEM_INLINE quatf GEM_VECTORCALL operator/(const quatf& lhs, const float rhs)
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

    GEM_INLINE quatf GEM_VECTORCALL operator*(const float lhs, const quatf& rhs)
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

#pragma region quatd

    struct quatd
    {
        union
        {
            struct
            {
                double x, y, z, w;
            };

            struct
            {
                double3 v;
                double s;
            };

            double components[4];
        };

        static quatd identity();

        static quatd GEM_VECTORCALL rotate_euler(const double pitch, const double yaw, const double roll);

        static quatd GEM_VECTORCALL rotate_axis_angle(const double3& axis, const double angle);

        static quatd GEM_VECTORCALL rotate_matrix(const double3x3& m);

        quatd(const double x, const double y, const double z, const double w);

        quatd(const double* o);

        quatd(const quatd& o);

        quatd() = default;

        double length_squared() const;

        double length() const;

        quatd& normalize();

        quatd& safe_normalize(const double tolerance = 0.001f);

        quatd& conjugate();

        quatd& inverse();

        double3x3 matrix3x3() const;

        double4x3 matrix4x3() const;

        double4x4 matrix4x4() const;

        double3 euler() const;

        double3 GEM_VECTORCALL transform_point(const double3& p) const;

        quatd& GEM_VECTORCALL operator=(const quatd& rhs);

        quatd& GEM_VECTORCALL operator+=(const quatd& rhs);

        quatd& GEM_VECTORCALL operator-=(const quatd& rhs);

        quatd& GEM_VECTORCALL operator*=(const quatd& rhs);

        quatd& GEM_VECTORCALL operator*=(const double rhs);

        quatd& GEM_VECTORCALL operator/=(const double rhs);

        quatd operator-() const;
    };

    double GEM_VECTORCALL dot(const quatd& lhs, const quatd& rhs);

    double GEM_VECTORCALL length(const quatd& rhs);

    double GEM_VECTORCALL length_squared(const quatd& val);

    quatd GEM_VECTORCALL normalize(const quatd& rhs);

    quatd GEM_VECTORCALL safe_normalize(const quatd& rhs, const double tolerance = 0.001f);

    quatd GEM_VECTORCALL conjugate(const quatd& rhs);

    quatd GEM_VECTORCALL inverse(const quatd& rhs);

    quatd GEM_VECTORCALL slerp(quatd q0, quatd q1, const double t);

    quatd GEM_VECTORCALL operator+(const quatd& lhs, const quatd& rhs);

    quatd GEM_VECTORCALL operator-(const quatd& lhs, const quatd& rhs);

    quatd GEM_VECTORCALL operator*(const quatd& lhs, const quatd& rhs);

    double3 GEM_VECTORCALL operator*(const double3& lhs, const quatd& rhs);

    double3 GEM_VECTORCALL operator*(const quatd& lhs, const double3& rhs);

    quatd GEM_VECTORCALL operator*(const quatd& lhs, const double rhs);

    quatd GEM_VECTORCALL operator*(const double lhs, const quatd& rhs);

    quatd GEM_VECTORCALL operator/(const quatd& lhs, const double rhs);

    GEM_INLINE quatd quatd::identity()
    {
        return { 0, 0, 0, 1 };
    }

    GEM_INLINE quatd GEM_VECTORCALL quatd::rotate_euler(const double pitch, const double yaw, const double roll)
    {
        double halfRoll = roll * 0.5f;
        double halfPitch = pitch * 0.5f;
        double halfYaw = yaw * 0.5f;
        double cosHalfRoll = cos(halfRoll);
        double cosHalfPitch = cos(halfPitch);
        double cosHalfYaw = cos(halfYaw);
        double sinHalfRoll = sin(halfRoll);
        double sinHalfPitch = sin(halfPitch);
        double sinHalfYaw = sin(halfYaw);

        return
        {
            (cosHalfYaw * sinHalfPitch * cosHalfRoll) + (sinHalfYaw * cosHalfPitch * sinHalfRoll),
            (sinHalfYaw * cosHalfPitch * cosHalfRoll) - (cosHalfYaw * sinHalfPitch * sinHalfRoll),
            (cosHalfYaw * cosHalfPitch * sinHalfRoll) - (sinHalfYaw * sinHalfPitch * cosHalfRoll),
            (cosHalfYaw * cosHalfPitch * cosHalfRoll) + (sinHalfYaw * sinHalfPitch * sinHalfRoll)
        };
    }

    GEM_INLINE quatd GEM_VECTORCALL quatd::rotate_axis_angle(const double3& axis, const double angle)
    {
        double a = angle * 0.5f;
        double s = sin(a);
        double c = cos(a);

        return
        {
            axis.x * s,
            axis.y * s,
            axis.z * s,
            c
        };
    }

    GEM_INLINE quatd GEM_VECTORCALL quatd::rotate_matrix(const double3x3& m)
    {
        quatd o;
        const double m11 = m.m00, m12 = m.m01, m13 = m.m02;
        const double m21 = m.m10, m22 = m.m11, m23 = m.m12;
        const double m31 = m.m20, m32 = m.m21, m33 = m.m22;

        // Determine which of w, x, y or z has the largest absolute value
        double fourWSquaredMinus1 = +m11 + m22 + m33;
        double fourXSquaredMinus1 = +m11 - m22 - m33;
        double fourYSquaredMinus1 = -m11 + m22 - m33;
        double fourZSquaredMinus1 = -m11 - m22 + m33;

        int biggestIndex = 0;
        double fourBiggestSquardeMinus1 = fourWSquaredMinus1;
        if (fourXSquaredMinus1 > fourBiggestSquardeMinus1)
        {
            fourBiggestSquardeMinus1 = fourXSquaredMinus1;
            biggestIndex = 1;
        }
        if (fourYSquaredMinus1 > fourBiggestSquardeMinus1)
        {
            fourBiggestSquardeMinus1 = fourYSquaredMinus1;
            biggestIndex = 2;
        }
        if (fourZSquaredMinus1 > fourBiggestSquardeMinus1)
        {
            fourBiggestSquardeMinus1 = fourZSquaredMinus1;
            biggestIndex = 3;
        }

        double biggestVal = sqrt(fourBiggestSquardeMinus1 + 1) * .5;
        double mult = 0.25 / biggestVal;

        switch (biggestIndex)
        {
        case 0:
        {
            o.x = (m23 - m32) * mult;
            o.y = (m31 - m13) * mult;
            o.z = (m12 - m21) * mult;
            o.w = biggestVal;
            break;
        }
        case 1:
        {
            o.x = biggestVal;
            o.y = (m12 + m21) * mult;
            o.z = (m31 + m13) * mult;
            o.w = (m23 - m32) * mult;
            break;
        }
        case 2:
        {
            o.x = (m12 + m21) * mult;
            o.y = biggestVal;
            o.z = (m23 + m32) * mult;
            o.w = (m31 - m13) * mult;
            break;
        }
        case 3:
        {
            o.x = (m31 + m13) * mult;
            o.y = (m23 + m32) * mult;
            o.z = biggestVal;
            o.w = (m12 - m21) * mult;
            break;
        }
        default:
            o.x = 0.0;
            o.y = 0.0;
            o.z = 0.0;
            o.w = 1.0;
            break;
        }
        return o;
    }

    GEM_INLINE quatd::quatd(const double x, const double y, const double z, const double w)
        : x(x)
        , y(y)
        , z(z)
        , w(w)
    {

    }

    GEM_INLINE quatd::quatd(const double* o)
        : x(o[0])
        , y(o[1])
        , z(o[2])
        , w(o[3])
    {

    }

    GEM_INLINE quatd::quatd(const quatd& o)
        : x(o.x)
        , y(o.y)
        , z(o.z)
        , w(o.w)
    {

    }

    GEM_INLINE double quatd::length_squared() const
    {
        return (x * x) + (y * y) + (z * z) + (w * w);
    }

    GEM_INLINE double quatd::length() const
    {
        return sqrt((x * x) + (y * y) + (z * z) + (w * w));
    }

    GEM_INLINE quatd& quatd::normalize()
    {
        double il = 1.0 / sqrt((x * x) + (y * y) + (z * z) + (w * w));
        x *= il;
        y *= il;
        z *= il;
        w *= il;
        return *this;
    }

    GEM_INLINE quatd& quatd::safe_normalize(const double tolerance /*= 0.001f*/)
    {
        double l = sqrt((x * x) + (y * y) + (z * z) + (w * w));
        if (l > tolerance)
        {
            double il = 1.0 / l;
            x *= il;
            y *= il;
            z *= il;
            w *= il;
        }

        return *this;
    }

    GEM_INLINE quatd& quatd::conjugate()
    {
        x = -x;
        y = -y;
        z = -z;
        return *this;
    }

    GEM_INLINE quatd& quatd::inverse()
    {
        double il = 1.0 / sqrt((x * x) + (y * y) + (z * z) + (w * w));
        x = -x * il;
        y = -y * il;
        z = -z * il;
        w = w * il;
        return *this;
    }

    GEM_INLINE double3x3 quatd::matrix3x3() const
    {
        double x2 = x * x;
        double y2 = y * y;
        double z2 = z * z;
        double wx = w * x;
        double wy = w * y;
        double wz = w * z;
        double xy = x * y;
        double xz = x * z;
        double yz = y * z;
        return
        {
            1.0 - (2.0 * y2) - (2.0 * z2),
            (2.0 * xy) + (2.0 * wz),
            (2.0 * xz) - (2.0 * wy),
            (2.0 * xy) - (2.0 * wz),
            1.0 - (2.0 * x2) - (2.0 * z2),
            (2.0 * yz) + (2.0 * wx),
            (2.0 * xz) + (2.0 * wy),
            (2.0 * yz) - (2.0 * wx),
            1.0 - (2.0 * x2) - (2.0 * y2)
        };
    }

    GEM_INLINE double4x3 quatd::matrix4x3() const
    {
        double x2 = x * x;
        double y2 = y * y;
        double z2 = z * z;
        double wx = w * x;
        double wy = w * y;
        double wz = w * z;
        double xy = x * y;
        double xz = x * z;
        double yz = y * z;

        return
        {
            1.0 - (2.0 * y2) - (2.0 * z2),
            (2.0 * xy) + (2.0 * wz),
            (2.0 * xz) - (2.0 * wy),
            (2.0 * xy) - (2.0 * wz),
            1.0 - (2.0 * x2) - (2.0 * z2),
            (2.0 * yz) + (2.0 * wx),
            (2.0 * xz) + (2.0 * wy),
            (2.0 * yz) - (2.0 * wx),
            1.0 - (2.0 * x2) - (2.0 * y2),
            0.0,
            0.0,
            0.0
        };
    }

    GEM_INLINE double4x4 quatd::matrix4x4() const
    {
        double x2 = x * x;
        double y2 = y * y;
        double z2 = z * z;
        double wx = w * x;
        double wy = w * y;
        double wz = w * z;
        double xy = x * y;
        double xz = x * z;
        double yz = y * z;

        return
        {
            1.0 - (2.0 * y2) - (2.0 * z2),
            (2.0 * xy) + (2.0 * wz),
            (2.0 * xz) - (2.0 * wy),
            0.0,
            (2.0 * xy) - (2.0 * wz),
            1.0 - (2.0 * x2) - (2.0 * z2),
            (2.0 * yz) + (2.0 * wx),
            0.0,
            (2.0 * xz) + (2.0 * wy),
            (2.0 * yz) - (2.0 * wx),
            1.0 - (2.0 * x2) - (2.0 * y2),
            0.0,
            0.0,
            0.0,
            0.0,
            1.0
        };
    }

    GEM_INLINE double3 quatd::euler() const
    {
        double3 o;
        auto sp = -2.0 * ((y * z) - (w * x));
        if (fabs(sp) > 0.9999)
        {
            o.x = 1.570796f * sp;
            o.y = atan2((-x * z) + (w * y), 0.5f - (y * y) - (z * z));
            o.z = 0.0;
        }
        else
        {
            o.x = asin(sp);
            o.y = atan2((x * z) + (w * y), 0.5f - (x * x) - (y * y));
            o.z = atan2((x * y) + (w * z), 0.5f - (x * x) - (z * z));
        }
        return o;
    }

    GEM_INLINE double3 GEM_VECTORCALL quatd::transform_point(const double3& p) const
    {
        // qvq* = (b + c)(v + 0)(-b + c)
        //      = v + 2c(bxv) + 2(bxbxv)
        double3 bxv = cross(v, p);
        double3 bxbxv = cross(v, bxv);
        return
        {
            p.x + ((bxv.x * w) + bxbxv.x) * 2.0,
            p.y + ((bxv.y * w) + bxbxv.y) * 2.0,
            p.z + ((bxv.z * w) + bxbxv.z) * 2.0
        };
    }

    GEM_INLINE quatd& GEM_VECTORCALL quatd::operator=(const quatd& rhs)
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        w = rhs.w;
        return *this;
    }

    GEM_INLINE quatd& GEM_VECTORCALL quatd::operator+=(const quatd& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        w += rhs.w;
        return *this;
    }

    GEM_INLINE quatd& GEM_VECTORCALL quatd::operator-=(const quatd& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        w -= rhs.w;
        return *this;
    }

    GEM_INLINE quatd& GEM_VECTORCALL quatd::operator*=(const quatd& rhs)
    {
        double ow = (w * rhs.w) - ((x * rhs.x) + (y * rhs.y) + (z * rhs.z));
        double ox = (w * rhs.x) + (rhs.w * x) + ((y * rhs.z) - (z * rhs.y));
        double oy = (w * rhs.y) + (rhs.w * y) + ((z * rhs.x) - (x * rhs.z));
        double oz = (w * rhs.z) + (rhs.w * z) + ((x * rhs.y) - (y * rhs.x));
        x = ox;
        y = oy;
        z = oz;
        w = ow;
        return *this;
    }

    GEM_INLINE quatd& GEM_VECTORCALL quatd::operator*=(const double rhs)
    {
        x *= rhs;
        y *= rhs;
        z *= rhs;
        w *= rhs;
        return *this;
    }

    GEM_INLINE quatd& GEM_VECTORCALL quatd::operator/=(const double rhs)
    {
        double inv = 1.0 / rhs;
        x *= inv;
        y *= inv;
        z *= inv;
        w *= inv;
        return *this;
    }

    GEM_INLINE quatd quatd::operator-() const
    {
        quatd out;
        out.x = -x;
        out.y = -y;
        out.z = -z;
        out.w = -w;
        return out;
    }

    GEM_INLINE double GEM_VECTORCALL dot(const quatd& lhs, const quatd& rhs)
    {
        return (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z) + (lhs.w * rhs.w);
    }

    GEM_INLINE double GEM_VECTORCALL length_squared(const quatd& val)
    {
        return (val.x * val.x) + (val.y * val.y) + (val.z * val.z) + (val.w * val.w);
    }

    GEM_INLINE double GEM_VECTORCALL length(const quatd& rhs)
    {
        return sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w));
    }

    GEM_INLINE quatd GEM_VECTORCALL normalize(const quatd& rhs)
    {
        double il = 1.0 / sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w));
        return
        {
            rhs.x * il,
            rhs.y * il,
            rhs.z * il,
            rhs.w * il
        };
    }

    GEM_INLINE quatd GEM_VECTORCALL safe_normalize(const quatd& rhs, const double tolerance /*= 0.001f*/)
    {
        double l = sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w));
        double il = (l > tolerance) ? (1.0 / l) : 1.0;
        return
        {
            rhs.x * il,
            rhs.y * il,
            rhs.z * il,
            rhs.w * il
        };
    }

    GEM_INLINE quatd GEM_VECTORCALL conjugate(const quatd& rhs)
    {
        return { -rhs.x, -rhs.y, -rhs.z, rhs.w };
    }

    GEM_INLINE quatd GEM_VECTORCALL inverse(const quatd& rhs)
    {
        double il = 1.0 / sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w));
        return { -rhs.x * il, -rhs.y * il, -rhs.z * il, rhs.w * il };
    }

    GEM_INLINE quatd GEM_VECTORCALL slerp(quatd q0, quatd q1, const double t)
    {
        double cosAngle = (q0.x * q1.x) + (q0.y * q1.y) + (q0.z * q1.z) + (q0.w * q1.w);
        if (cosAngle < 0.0) {
            q1 = -q1;
            cosAngle = -cosAngle;
        }

        double k0, k1;

        // Check for divide by zero
        if (cosAngle > 0.9999) {
            k0 = 1.0 - t;
            k1 = t;
        }
        else {
            double angle = acos(cosAngle);
            double oneOverSinAngle = 1.0 / sin(angle);

            k0 = ((sin(1.0 - t) * angle) * oneOverSinAngle);
            k1 = (sin(t * angle) * oneOverSinAngle);
        }

        q0 = q0 * k0;
        q1 = q1 * k1;

        return q0 + q1;
    }

    GEM_INLINE quatd GEM_VECTORCALL operator+(const quatd& lhs, const quatd& rhs)
    {
        return
        {
            lhs.x + rhs.x,
            lhs.y + rhs.y,
            lhs.z + rhs.z,
            lhs.w + rhs.w
        };
    }

    GEM_INLINE quatd GEM_VECTORCALL operator-(const quatd& lhs, const quatd& rhs)
    {
        return
        {
            lhs.x - rhs.x,
            lhs.y - rhs.y,
            lhs.z - rhs.z,
            lhs.w - rhs.w
        };
    }

    GEM_INLINE quatd GEM_VECTORCALL operator*(const quatd& lhs, const quatd& rhs)
    {
        return
        {
            (lhs.w * rhs.x) + (lhs.x * rhs.w) + ((lhs.y * rhs.z) - (lhs.z * rhs.y)),
            (lhs.w * rhs.y) + (lhs.y * rhs.w) + ((lhs.z * rhs.x) - (lhs.x * rhs.z)),
            (lhs.w * rhs.z) + (lhs.z * rhs.w) + ((lhs.x * rhs.y) - (lhs.y * rhs.x)),
            (lhs.w * rhs.w) - ((lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z))
        };
    }

    GEM_INLINE quatd GEM_VECTORCALL operator*(const quatd& lhs, const double rhs)
    {
        return
        {
            lhs.x * rhs,
            lhs.y * rhs,
            lhs.z * rhs,
            lhs.w * rhs
        };
    }

    GEM_INLINE quatd GEM_VECTORCALL operator/(const quatd& lhs, const double rhs)
    {
        double inv = 1.0 / rhs;
        return
        {
            lhs.x * inv,
            lhs.y * inv,
            lhs.z * inv,
            lhs.w * inv
        };
    }

    GEM_INLINE quatd GEM_VECTORCALL operator*(const double lhs, const quatd& rhs)
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