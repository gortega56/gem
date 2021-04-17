#pragma once
#include "common/defines.h"
#include "matrix.h"

namespace gem
{
#pragma region quatf

    struct quatf
    {
        union
        {
            struct
            {
                float x, y, z, w;
            };

            struct
            {
                float3 v;
                float s;
            };

            float components[4];
        };

        static quatf identity();

        static quatf GEM_VECTORCALL rotate_euler(const float pitch, const float yaw, const float roll);

        static quatf GEM_VECTORCALL rotate_axis_angle(const float3& axis, const float angle);

        static quatf GEM_VECTORCALL rotate_matrix(const float3x3& m);

        quatf(const float x, const float y, const float z, const float w);

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

    float3 GEM_VECTORCALL operator*(const float3& lhs, const quatf& rhs);

    float3 GEM_VECTORCALL operator*(const quatf& lhs, const float3& rhs);

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
        const float m11 = m.entries1d[0], m12 = m.entries1d[1], m13 = m.entries1d[2];
        const float m21 = m.entries1d[3], m22 = m.entries1d[4], m23 = m.entries1d[5];
        const float m31 = m.entries1d[6], m32 = m.entries1d[7], m33 = m.entries1d[8];

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

    GEM_INLINE quatf::quatf(const float x, const float y, const float z, const float w)
        : x(x)
        , y(y)
        , z(z)
        , w(w)
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
        quatf out;
        out.x = -x;
        out.y = -y;
        out.z = -z;
        out.w = -w;
        return out;
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
        return sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w));
    }

    GEM_INLINE quatf GEM_VECTORCALL normalize(const quatf& rhs)
    {
        float il = 1.0f / sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w));
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

    GEM_INLINE float3 GEM_VECTORCALL operator*(const float3& lhs, const quatf& rhs)
    {
        return rhs * lhs;
    }

    GEM_INLINE float3 GEM_VECTORCALL operator*(const quatf& lhs, const float3& rhs)
    {
        float3 vxp = cross(lhs.v, rhs);
        float3 vxpxv = cross(lhs.v, vxp);
        return
        {
            rhs.x + ((vxp.x * lhs.w) + vxpxv.x) * 2.0f,
            rhs.y + ((vxp.y * lhs.w) + vxpxv.y) * 2.0f,
            rhs.z + ((vxp.z * lhs.w) + vxpxv.z) * 2.0f
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
        const double m11 = m.entries1d[0], m12 = m.entries1d[1], m13 = m.entries1d[2];
        const double m21 = m.entries1d[3], m22 = m.entries1d[4], m23 = m.entries1d[5];
        const double m31 = m.entries1d[6], m32 = m.entries1d[7], m33 = m.entries1d[8];

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

        double biggestVal = sqrt(fourBiggestSquardeMinus1 + 1) * .5f;
        double mult = 0.25f / biggestVal;

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
        double il = 1.0f / sqrt((x * x) + (y * y) + (z * z) + (w * w));
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
            double il = 1.0f / l;
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
        double il = 1.0f / sqrt((x * x) + (y * y) + (z * z) + (w * w));
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

    GEM_INLINE double3 quatd::euler() const
    {
        double3 o;
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
        double inv = 1.0f / rhs;
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
        double il = 1.0f / sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w));
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
        double il = (l > tolerance) ? (1.0f / l) : 1.0f;
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
        double il = 1.0f / sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w));
        return { -rhs.x * il, -rhs.y * il, -rhs.z * il, rhs.w * il };
    }

    GEM_INLINE quatd GEM_VECTORCALL slerp(quatd q0, quatd q1, const double t)
    {
        double cosAngle = (q0.x * q1.x) + (q0.y * q1.y) + (q0.z * q1.z) + (q0.w * q1.w);
        if (cosAngle < 0.0f) {
            q1 = -q1;
            cosAngle = -cosAngle;
        }

        double k0, k1;

        // Check for divide by zero
        if (cosAngle > 0.9999f) {
            k0 = 1.0f - t;
            k1 = t;
        }
        else {
            double angle = acos(cosAngle);
            double oneOverSinAngle = 1.0f / sin(angle);

            k0 = ((sin(1.0f - t) * angle) * oneOverSinAngle);
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

    GEM_INLINE double3 GEM_VECTORCALL operator*(const double3& lhs, const quatd& rhs)
    {
        return rhs * lhs;
    }

    GEM_INLINE double3 GEM_VECTORCALL operator*(const quatd& lhs, const double3& rhs)
    {
        double3 vxp = cross(lhs.v, rhs);
        double3 vxpxv = cross(lhs.v, vxp);
        return
        {
            rhs.x + ((vxp.x * lhs.w) + vxpxv.x) * 2.0f,
            rhs.y + ((vxp.y * lhs.w) + vxpxv.y) * 2.0f,
            rhs.z + ((vxp.z * lhs.w) + vxpxv.z) * 2.0f
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
        double inv = 1.0f / rhs;
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