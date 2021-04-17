#pragma once
#include "common/defines.h"
#include "vector.h"

namespace gem
{
#pragma region float2x2

    struct float2x2
    {
        union
        {
            struct
            {
                float m00, m01;
                float m10, m11;
            };

            struct
            {
                float2 u, v;
            };

            float2 rows[2];
            float entries1d[4];
            float entries2d[2][2];
        };

        float2x2(const float u1, const float u2, const float v1, const float v2);

        float2x2(const float* p_entries);

        float2x2(const float2& u, const float2& v);

        float2x2(const float2x2& other);

        float2x2() = default;

        float2x2 transpose() const;

        float2x2 inverse() const;

        float determinant()	const;

        float2x2& GEM_VECTORCALL operator+=(const float2x2& rhs);

        float2x2& GEM_VECTORCALL operator-=(const float2x2& rhs);

        float2x2& GEM_VECTORCALL operator*=(const float rhs);

        float2x2& GEM_VECTORCALL operator=(const float2x2& rhs);

        float2& GEM_VECTORCALL operator[](const unsigned int index);

        float2 GEM_VECTORCALL operator[](const unsigned int index) const;

        float GEM_VECTORCALL operator()(const unsigned int row, const unsigned int column);

        float2x2 operator-() const;
    };

    float2x2 GEM_VECTORCALL operator+(const float2x2& lhs, const float2x2& rhs);

    float2x2 GEM_VECTORCALL operator-(const float2x2& lhs, const float2x2& rhs);

    float2x2 GEM_VECTORCALL operator*(const float2x2& lhs, const float2x2& rhs);

    float2x2 GEM_VECTORCALL operator*(const float2x2& lhs, const float rhs);

    float2x2 GEM_VECTORCALL operator*(const float lhs, const float2x2& rhs);

    float2 GEM_VECTORCALL operator*(const float2& lhs, const float2x2& rhs);

    GEM_INLINE float2x2::float2x2(const float u1, const float u2, const float v1, const float v2)
        : m00(u1), m01(u2)
        , m10(v1), m11(v2)
    {

    }

    GEM_INLINE float2x2::float2x2(const float* p_entries)
        : m00(p_entries[0]), m01(p_entries[1])
        , m10(p_entries[2]), m11(p_entries[3])
    {

    }

    GEM_INLINE float2x2::float2x2(const float2& u, const float2& v)
        : m00(u.x), m01(u.y)
        , m10(v.x), m11(v.y)
    {

    }

    GEM_INLINE float2x2::float2x2(const float2x2& o)
        : m00(o.m00), m01(o.m01)
        , m10(o.m10), m11(o.m11)
    {

    }

    GEM_INLINE float2x2 float2x2::transpose() const
    {
        return
        {
            m00, m10,
            m01, m11
        };
    }

    GEM_INLINE float2x2 float2x2::inverse() const
    {
        float2x2 out;
        float det = (entries1d[0] * entries1d[3]) - (entries1d[1] * entries1d[2]);
        float inv = 1.0f / det;
        out.entries1d[0] = +entries1d[3] * inv;
        out.entries1d[1] = -entries1d[1] * inv;
        out.entries1d[2] = -entries1d[2] * inv;
        out.entries1d[3] = +entries1d[0] * inv;
        return out;
    }

    GEM_INLINE float float2x2::determinant() const
    {
        //     0 1      
        // det 2 3 =    (0 * 3) - (1 - 2)     
        return (entries1d[0] * entries1d[3]) - (entries1d[1] * entries1d[2]);
    }

    GEM_INLINE float2x2& GEM_VECTORCALL float2x2::operator+=(const float2x2& rhs)
    {
        m00 += rhs.m00; m01 += rhs.m01;
        m10 += rhs.m10; m11 += rhs.m11;
        return *this;
    }

    GEM_INLINE float2x2& GEM_VECTORCALL float2x2::operator-=(const float2x2& rhs)
    {
        m00 -= rhs.m00; m01 -= rhs.m01;
        m10 -= rhs.m10; m11 -= rhs.m11;
        return *this;
    }

    GEM_INLINE float2x2& GEM_VECTORCALL float2x2::operator*=(const float rhs)
    {
        m00 *= rhs; m01 *= rhs;
        m10 *= rhs; m11 *= rhs;
        return *this;
    }

    GEM_INLINE float2x2& GEM_VECTORCALL float2x2::operator=(const float2x2& rhs)
    {
        m00 = rhs.m00; m01 = rhs.m01;
        m10 = rhs.m10; m11 = rhs.m11;
        return *this;
    }

    GEM_INLINE float2x2 float2x2::operator-() const
    {
        return
        {
            -m00, -m01,
            -m10, -m11
        };
    }

    GEM_INLINE float2& GEM_VECTORCALL float2x2::operator[](const unsigned int index)
    {
        return rows[index];
    }

    GEM_INLINE float2 GEM_VECTORCALL float2x2::operator[](const unsigned int index) const
    {
        return rows[index];
    }

    GEM_INLINE float GEM_VECTORCALL float2x2::operator()(const unsigned int row, const unsigned int column)
    {
        return entries2d[row][column];
    }

    GEM_INLINE float2x2 GEM_VECTORCALL operator+(const float2x2& lhs, const float2x2& rhs)
    {
        return
        {
            lhs.entries1d[0] + rhs.entries1d[0], lhs.entries1d[1] + rhs.entries1d[1],
            lhs.entries1d[2] + rhs.entries1d[2], lhs.entries1d[3] + rhs.entries1d[3]
        };
    }

    GEM_INLINE float2x2 GEM_VECTORCALL operator-(const float2x2& lhs, const float2x2& rhs)
    {
        return
        {
            lhs.entries1d[0] - rhs.entries1d[0], lhs.entries1d[1] - rhs.entries1d[1],
            lhs.entries1d[2] - rhs.entries1d[2], lhs.entries1d[3] - rhs.entries1d[3]
        };
    }

    GEM_INLINE float2x2 GEM_VECTORCALL operator*(const float2x2& lhs, const float2x2& rhs)
    {
        return
        {
            (lhs.entries1d[0] * rhs.entries1d[0]) + (lhs.entries1d[1] * rhs.entries1d[2]),
            (lhs.entries1d[0] * rhs.entries1d[1]) + (lhs.entries1d[1] * rhs.entries1d[3]),
            (lhs.entries1d[2] * rhs.entries1d[0]) + (lhs.entries1d[3] * rhs.entries1d[2]),
            (lhs.entries1d[2] * rhs.entries1d[1]) + (lhs.entries1d[3] * rhs.entries1d[3])
        };
    }

    GEM_INLINE float2 GEM_VECTORCALL operator*(const float2& lhs, const float2x2& rhs)
    {
        return
        {
            (lhs.components[0] * rhs.entries1d[0]) + (lhs.components[1] * rhs.entries1d[2]),
            (lhs.components[0] * rhs.entries1d[1]) + (lhs.components[1] * rhs.entries1d[3])
        };
    }

    GEM_INLINE float2x2 GEM_VECTORCALL operator*(const float lhs, const float2x2& rhs)
    {
        return
        {
            rhs.entries1d[0] * lhs, rhs.entries1d[1] * lhs,
            rhs.entries1d[2] * lhs, rhs.entries1d[3] * lhs
        };
    }

    GEM_INLINE float2x2 GEM_VECTORCALL operator*(const float2x2& lhs, const float rhs)
    {
        return
        {
            lhs.entries1d[0] * rhs, lhs.entries1d[1] * rhs,
            lhs.entries1d[2] * rhs, lhs.entries1d[3] * rhs
        };
    }

#pragma endregion

#pragma region float3x3

    struct float3x3
    {
        union
        {
            struct
            {
                float m00, m01, m02;
                float m10, m11, m12;
                float m20, m21, m22;
            };

            struct
            {
                float3 u, v, w;
            };

            float3 rows[3];
            float entries1d[9];
            float entries2d[3][3];
        };

        static float3x3 identity();

        static float3x3 GEM_VECTORCALL rotate_euler(const float pitch, const float yaw, const float roll);

        static float3x3 GEM_VECTORCALL rotate_axis_angle(const float3& a, const float radians);

        static float3x3 GEM_VECTORCALL rotate_axis_x(const float radians);

        static float3x3 GEM_VECTORCALL rotate_axis_y(const float radians);

        static float3x3 GEM_VECTORCALL rotate_axis_z(const float radians);

        static float3x3 GEM_VECTORCALL rotate_forward_lh(const float3& direction, const float3& up);

        static float3x3 GEM_VECTORCALL rotate_forward_rh(const float3& direction, const float3& up);

        static float3x3 GEM_VECTORCALL scale(const float3& s);

        float3x3(const float u1, const float u2, const float u3, const float v1, const float v2, const float v3, const float w1, const float w2, const float w3);

        float3x3(const float* p_entries);

        float3x3(const float3& u, const float3& v, const float3& w);

        float3x3(const float3x3& o);

        float3x3() = default;

        float3x3 transpose() const;

        float3x3 inverse() const;

        float determinant()	const;

        float3 euler() const;

        float3x3& GEM_VECTORCALL operator=(const float3x3& rhs);

        float3x3& GEM_VECTORCALL operator+=(const float3x3& rhs);

        float3x3& GEM_VECTORCALL operator-=(const float3x3& rhs);

        float3x3& GEM_VECTORCALL operator*=(const float rhs);

        float3& GEM_VECTORCALL operator[](const unsigned int index);

        float3 GEM_VECTORCALL operator[](const unsigned int index) const;

        float GEM_VECTORCALL operator()(const unsigned int row, const unsigned int column);
    
        float3x3 operator-() const;
    };

    float3x3 GEM_VECTORCALL operator+(const float3x3& lhs, const float3x3& rhs);

    float3x3 GEM_VECTORCALL operator-(const float3x3& lhs, const float3x3& rhs);

    float3x3 GEM_VECTORCALL operator*(const float3x3& lhs, const float3x3& rhs);

    float3x3 GEM_VECTORCALL operator*(const float3x3& lhs, const float rhs);

    float3x3 GEM_VECTORCALL operator*(const float lhs, const float3x3& rhs);

    float3 GEM_VECTORCALL operator*(const float3& lhs, const float3x3& rhs);

    GEM_INLINE float3x3 float3x3::identity()
    {
        return
        {
            1.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 1.0f
        };
    }

    GEM_INLINE float3x3 GEM_VECTORCALL float3x3::rotate_euler(const float pitch, const float yaw, const float roll)
    {
        float ch = cosf(yaw); float cp = cosf(pitch); float cb = cosf(roll);
        float sh = sinf(yaw); float sp = sinf(pitch); float sb = sinf(roll);

        return
        {
             ch * cb + sh * sp * sb, sb * cp, -sh * cb + ch * sp * sb,
            -ch * sb + sh * sp * cb, cb * cp,  sb * sh + ch * sp * cb,
             sh * cp, -sp,  ch * cp
        };
    }

    GEM_INLINE float3x3 GEM_VECTORCALL float3x3::rotate_axis_angle(const float3& a, const float radians)
    {
        float s = sinf(radians);
        float c = cosf(radians);

        return
        {
            c + ((1 - c) * a.x * a.x), ((1 - c) * a.x * a.y) + (s * a.z), ((1 - c) * a.x * a.z) - (s * a.y),
            ((1 - c) * a.x * a.y) - (s * a.z), c + ((1 - c) * a.y * a.y), ((1 - c) * a.y * a.z) + (s * a.x),
            ((1 - c) * a.x * a.z) + (s * a.y), ((1 - c) * a.y * a.z) - (s * a.x), c + ((1 - c) * a.z * a.z)
        };
    }

    GEM_INLINE float3x3 GEM_VECTORCALL float3x3::rotate_axis_x(const float radians)
    {
        float s = sinf(radians);
        float c = cosf(radians);
        return
        {
            1.0f, 0.0f, 0.0f,
            0.0f, c, s,
            0.0f, -s, c
        };
    }

    GEM_INLINE float3x3 GEM_VECTORCALL float3x3::rotate_axis_y(const float radians)
    {
        float s = sin(radians);
        float c = cos(radians);
        return
        {
            c, 0.0f, -s,
            0.0f, 1.0f, 0.0f,
            s, 0.0f, c
        };
    }

    GEM_INLINE float3x3 GEM_VECTORCALL float3x3::rotate_axis_z(const float radians)
    {
        float s = sin(radians);
        float c = cos(radians);
        return
        {
            c, s, 0.0f,
            -s, c, 0.0f,
            0.0f, 0.0f, 1.0f
        };
    }

    GEM_INLINE float3x3 GEM_VECTORCALL float3x3::rotate_forward_lh(const float3& direction, const float3& up)
    {
        float3 f = normalize(direction);
        float3 s = normalize(cross(up, f));
        float3 u = cross(f, s);

        return
        {
            s.x, s.y, s.z,
            u.x, u.y, u.z,
            f.x, f.y, f.z
        };
    }

    GEM_INLINE float3x3 GEM_VECTORCALL float3x3::rotate_forward_rh(const float3& direction, const float3& up)
    {
        return rotate_forward_lh(-direction, up);
    }

    GEM_INLINE float3x3 GEM_VECTORCALL float3x3::scale(const float3& s)
    {
        return
        {
            s.x, 0.0f, 0.0f,
            0.0f, s.y, 0.0f,
            0.0f, 0.0f, s.z
        };
    }

    GEM_INLINE float3x3::float3x3(
        const float u1, const float u2, const float u3,
        const float v1, const float v2, const float v3,
        const float w1, const float w2, const float w3)
        : m00(u1), m01(u2), m02(u3)
        , m10(v1), m11(v2), m12(v3)
        , m20(w1), m21(w2), m22(w3)
    {

    }

    GEM_INLINE float3x3::float3x3(const float* o)
        : m00(o[0]), m01(o[1]), m02(o[2])
        , m10(o[3]), m11(o[4]), m12(o[5])
        , m20(o[6]), m21(o[7]), m22(o[8])
    {

    }

    GEM_INLINE float3x3::float3x3(const float3& u, const float3& v, const float3& w)
        : m00(u.x), m01(u.y), m02(u.z)
        , m10(v.x), m11(v.y), m12(v.z)
        , m20(w.x), m21(w.y), m22(w.z)
    {

    };

    GEM_INLINE float3x3::float3x3(const float3x3& o)
        : m00(o.m00), m01(o.m01), m02(o.m02)
        , m10(o.m10), m11(o.m11), m12(o.m12)
        , m20(o.m20), m21(o.m21), m22(o.m22)
    {

    }

    GEM_INLINE float3x3 float3x3::transpose() const
    {
        return
        {
            m00, m10, m20,
            m01, m11, m21,
            m02, m12, m22
        };
    }

    GEM_INLINE float3x3 float3x3::inverse() const
    {
        float3x3 out;

        float det =
            (((entries1d[1] * entries1d[5]) - (entries1d[2] * entries1d[4])) * entries1d[6]) +
            (((entries1d[2] * entries1d[3]) - (entries1d[0] * entries1d[5])) * entries1d[7]) +
            (((entries1d[0] * entries1d[4]) - (entries1d[1] * entries1d[3])) * entries1d[8]);

        float c11 = +((entries1d[4] * entries1d[8]) - (entries1d[5] * entries1d[7]));
        float c12 = -((entries1d[3] * entries1d[8]) - (entries1d[5] * entries1d[6]));
        float c13 = +((entries1d[3] * entries1d[7]) - (entries1d[4] * entries1d[6]));
        float c21 = -((entries1d[1] * entries1d[8]) - (entries1d[2] * entries1d[7]));
        float c22 = +((entries1d[0] * entries1d[8]) - (entries1d[2] * entries1d[6]));
        float c23 = -((entries1d[0] * entries1d[7]) - (entries1d[1] * entries1d[6]));
        float c31 = +((entries1d[1] * entries1d[5]) - (entries1d[2] * entries1d[4]));
        float c32 = -((entries1d[0] * entries1d[5]) - (entries1d[2] * entries1d[3]));
        float c33 = +((entries1d[0] * entries1d[4]) - (entries1d[1] * entries1d[3]));

        float inv = 1.0f / det;
        out.entries1d[0] = c11 * inv;
        out.entries1d[1] = c21 * inv;
        out.entries1d[2] = c31 * inv;
        out.entries1d[3] = c12 * inv;
        out.entries1d[4] = c22 * inv;
        out.entries1d[5] = c32 * inv;
        out.entries1d[6] = c13 * inv;
        out.entries1d[7] = c23 * inv;
        out.entries1d[8] = c33 * inv;
        return out;
    }

    GEM_INLINE float float3x3::determinant() const
    {
        // ((row0 x row1) * row2)
        // (ay*bz - az*by) * cx
        // (az*bx - ax*bz) * cy
        // (ax*by - ay*bx) * cz

        return
            (((entries1d[1] * entries1d[5]) - (entries1d[2] * entries1d[4])) * entries1d[6]) +
            (((entries1d[2] * entries1d[3]) - (entries1d[0] * entries1d[5])) * entries1d[7]) +
            (((entries1d[0] * entries1d[4]) - (entries1d[1] * entries1d[3])) * entries1d[8]);
    }

    GEM_INLINE float3 float3x3::euler() const
    {
        float3 out;
        float sp = -entries1d[7];
        if (sp <= 1.0f)
        {
            out.x = -1.570796f;
        }
        else if (sp >= 1.0f)
        {
            out.x = -1.570796f;
        }
        else
        {
            out.x = asin(sp);
        }

        if (abs(sp) > 0.9999f)
        {
            out.z = 0.0f;
            out.y = atan2(-entries1d[2], entries1d[0]);
        }
        else
        {
            out.y = atan2(entries1d[6], entries1d[8]);
            out.z = atan2(entries1d[1], entries1d[4]);
        }
        return out;
    }

    GEM_INLINE float3x3& GEM_VECTORCALL float3x3::operator=(const float3x3& rhs)
    {
        m00 = rhs.m00; m01 = rhs.m01; m02 = rhs.m02; 
        m10 = rhs.m10; m11 = rhs.m11; m12 = rhs.m12; 
        m20 = rhs.m20; m21 = rhs.m21; m22 = rhs.m22;
        return *this;
    }

    GEM_INLINE float3x3& GEM_VECTORCALL float3x3::operator+=(const float3x3& rhs)
    {
        m00 += rhs.m00; m01 += rhs.m01; m02 += rhs.m02;
        m10 += rhs.m10; m11 += rhs.m11; m12 += rhs.m12;
        m20 += rhs.m20; m21 += rhs.m21; m22 += rhs.m22;
        return *this;
    }

    GEM_INLINE float3x3& GEM_VECTORCALL float3x3::operator-=(const float3x3& rhs)
    {
        m00 -= rhs.m00; m01 -= rhs.m01; m02 -= rhs.m02;
        m10 -= rhs.m10; m11 -= rhs.m11; m12 -= rhs.m12;
        m20 -= rhs.m20; m21 -= rhs.m21; m22 -= rhs.m22;
        return *this;
    }

    GEM_INLINE float3x3& GEM_VECTORCALL float3x3::operator*=(const float rhs)
    {
        m00 *= rhs; m01 *= rhs; m02 *= rhs;
        m10 *= rhs; m11 *= rhs; m12 *= rhs;
        m20 *= rhs; m21 *= rhs; m22 *= rhs;
        return *this;
    }

    GEM_INLINE float3& GEM_VECTORCALL float3x3::operator[](const unsigned int index)
    {
        return rows[index];
    }

    GEM_INLINE float3 GEM_VECTORCALL float3x3::operator[](const unsigned int index) const
    {
        return rows[index];
    }

    GEM_INLINE float GEM_VECTORCALL float3x3::operator()(const unsigned int row, const unsigned int column)
    {
        return entries2d[row][column];
    }

    GEM_INLINE float3x3 float3x3::operator-() const
    {
        return
        {
            -m00, -m01, -m02,
            -m10, -m11, -m12,
            -m20, -m21, -m22
        };
    }

    GEM_INLINE float3x3 GEM_VECTORCALL operator+(const float3x3& lhs, const float3x3& rhs)
    {
        return
        {
            lhs.entries1d[0] + rhs.entries1d[0], lhs.entries1d[1] + rhs.entries1d[1], lhs.entries1d[2] + rhs.entries1d[2],
            lhs.entries1d[3] + rhs.entries1d[3], lhs.entries1d[4] + rhs.entries1d[4], lhs.entries1d[5] + rhs.entries1d[5],
            lhs.entries1d[6] + rhs.entries1d[6], lhs.entries1d[7] + rhs.entries1d[7], lhs.entries1d[8] + rhs.entries1d[8],
        };
    }

    GEM_INLINE float3x3 GEM_VECTORCALL operator-(const float3x3& lhs, const float3x3& rhs)
    {
        return
        {
            lhs.entries1d[0] - rhs.entries1d[0], lhs.entries1d[1] - rhs.entries1d[1], lhs.entries1d[2] - rhs.entries1d[2],
            lhs.entries1d[3] - rhs.entries1d[3], lhs.entries1d[4] - rhs.entries1d[4], lhs.entries1d[5] - rhs.entries1d[5],
            lhs.entries1d[6] - rhs.entries1d[6], lhs.entries1d[7] - rhs.entries1d[7], lhs.entries1d[8] - rhs.entries1d[8],
        };
    }

    GEM_INLINE float3x3 GEM_VECTORCALL operator*(const float3x3& lhs, const float3x3& rhs)
    {
        return
        {
            (lhs.entries1d[0] * rhs.entries1d[0]) + (lhs.entries1d[1] * rhs.entries1d[3]) + (lhs.entries1d[2] * rhs.entries1d[6]),
            (lhs.entries1d[0] * rhs.entries1d[1]) + (lhs.entries1d[1] * rhs.entries1d[4]) + (lhs.entries1d[2] * rhs.entries1d[7]),
            (lhs.entries1d[0] * rhs.entries1d[2]) + (lhs.entries1d[1] * rhs.entries1d[5]) + (lhs.entries1d[2] * rhs.entries1d[8]),
            (lhs.entries1d[3] * rhs.entries1d[0]) + (lhs.entries1d[4] * rhs.entries1d[3]) + (lhs.entries1d[5] * rhs.entries1d[6]),
            (lhs.entries1d[3] * rhs.entries1d[1]) + (lhs.entries1d[4] * rhs.entries1d[4]) + (lhs.entries1d[5] * rhs.entries1d[7]),
            (lhs.entries1d[3] * rhs.entries1d[2]) + (lhs.entries1d[4] * rhs.entries1d[5]) + (lhs.entries1d[5] * rhs.entries1d[8]),
            (lhs.entries1d[6] * rhs.entries1d[0]) + (lhs.entries1d[7] * rhs.entries1d[3]) + (lhs.entries1d[8] * rhs.entries1d[6]),
            (lhs.entries1d[6] * rhs.entries1d[1]) + (lhs.entries1d[7] * rhs.entries1d[4]) + (lhs.entries1d[8] * rhs.entries1d[7]),
            (lhs.entries1d[6] * rhs.entries1d[2]) + (lhs.entries1d[7] * rhs.entries1d[5]) + (lhs.entries1d[8] * rhs.entries1d[8])
        };
    }

    GEM_INLINE float3 GEM_VECTORCALL operator*(const float3& lhs, const float3x3& rhs)
    {
        return
        {

            (lhs.components[0] * rhs.entries1d[0]) + (lhs.components[1] * rhs.entries1d[3]) + (lhs.components[2] * rhs.entries1d[6]),
            (lhs.components[0] * rhs.entries1d[1]) + (lhs.components[1] * rhs.entries1d[4]) + (lhs.components[2] * rhs.entries1d[7]),
            (lhs.components[0] * rhs.entries1d[2]) + (lhs.components[1] * rhs.entries1d[5]) + (lhs.components[2] * rhs.entries1d[8])
        };
    }

    GEM_INLINE float3x3 GEM_VECTORCALL operator*(const float lhs, const float3x3& rhs)
    {
        return
        {
            rhs.entries1d[0] * lhs, rhs.entries1d[1] * lhs, rhs.entries1d[2] * lhs,
            rhs.entries1d[3] * lhs, rhs.entries1d[4] * lhs, rhs.entries1d[5] * lhs,
            rhs.entries1d[6] * lhs, rhs.entries1d[7] * lhs, rhs.entries1d[8] * lhs
        };
    }

    GEM_INLINE float3x3 GEM_VECTORCALL operator*(const float3x3& lhs, const float rhs)
    {
        return
        {
            lhs.entries1d[0] * rhs, lhs.entries1d[1] * rhs, lhs.entries1d[2] * rhs,
            lhs.entries1d[3] * rhs, lhs.entries1d[4] * rhs, lhs.entries1d[5] * rhs,
            lhs.entries1d[6] * rhs, lhs.entries1d[7] * rhs, lhs.entries1d[8] * rhs
        };
    }

#pragma endregion

#pragma region float4x3

    struct float4x3
    {
        union
        {
            struct
            {
                float m00, m01, m02;
                float m10, m11, m12;
                float m20, m21, m22;
                float m30, m31, m32;
            };

            struct
            {
                float3 u, v, w, t;
            };

            float3 rows[4];
            float entries1d[12];
            float entries2d[4][3];
        };

        static float4x3 identity();

        static float4x3 GEM_VECTORCALL rotate_euler(const float pitch, const float yaw, const float roll);

        static float4x3 GEM_VECTORCALL rotate_axis_angle(const float3& a, const float angle);

        static float4x3 GEM_VECTORCALL rotate_axis_x(const float angle);

        static float4x3 GEM_VECTORCALL rotate_axis_y(const float angle);

        static float4x3 GEM_VECTORCALL rotate_axis_z(const float angle);

        static float4x3 GEM_VECTORCALL rotate_forward_lh(const float3& direction, const float3& position, const float3& up);

        static float4x3 GEM_VECTORCALL rotate_forward_rh(const float3& direction, const float3& position, const float3& up);

        static float4x3 GEM_VECTORCALL scale(const float3& s);

        static float4x3 GEM_VECTORCALL translate(const float3& t);

        float4x3(const float u1, const float u2, const float u3,
            const float v1, const float v2, const float v3,
            const float w1, const float w2, const float w3,
            const float t1, const float t2, const float t3);

        float4x3(const float* o);

        float4x3(const float3& u,
            const float3& v,
            const float3& w,
            const float3& t);

        float4x3(const float4x3& o);

        float4x3() = default;

        float4x3 transpose() const;

        float4x3 inverse() const;

        float determinant()	const;

        float3 euler() const;

        float4x3& GEM_VECTORCALL operator=(const float4x3& rhs);

        float4x3 GEM_VECTORCALL operator+=(const float4x3& rhs);

        float4x3 GEM_VECTORCALL operator-=(const float4x3& rhs);

        float4x3 GEM_VECTORCALL operator*=(const float rhs);

        float3& GEM_VECTORCALL operator[](const unsigned int index);

        float3 GEM_VECTORCALL operator[](const unsigned int index) const;

        float GEM_VECTORCALL operator()(const unsigned int row, const unsigned int column);

        float4x3 operator-() const;
    };

    float4x3 GEM_VECTORCALL operator+(const float4x3& lhs, const float4x3& rhs);

    float4x3 GEM_VECTORCALL operator-(const float4x3& lhs, const float4x3& rhs);

    float4x3 GEM_VECTORCALL operator*(const float4x3& lhs, const float4x3& rhs);

    float4x3 GEM_VECTORCALL operator*(const float4x3& lhs, const float rhs);

    float4x3 GEM_VECTORCALL operator*(const float lhs, const float4x3& rhs);

    float3 GEM_VECTORCALL operator*(const float3& lhs, const float4x3& rhs);

    GEM_INLINE float4x3 float4x3::identity()
    {
        return
        {
            1.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 1.0f,
            0.0f, 0.0f, 0.0f
        };
    }

    GEM_INLINE float4x3 GEM_VECTORCALL float4x3::rotate_euler(const float pitch, const float yaw, const float roll)
    {
        float ch = cosf(yaw); float cp = cosf(pitch); float cb = cosf(roll);
        float sh = sinf(yaw); float sp = sinf(pitch); float sb = sinf(roll);

        return
        {
             ch * cb + sh * sp * sb, sb * cp, -sh * cb + ch * sp * sb,
            -ch * sb + sh * sp * cb, cb * cp,  sb * sh + ch * sp * cb,
             sh * cp, -sp,  ch * cp,
             0.0f, 0.0f, 0.0f
        };
    }

    GEM_INLINE float4x3 GEM_VECTORCALL float4x3::rotate_axis_angle(const float3& a, const float angle)
    {
        float s = sin(angle);
        float c = cos(angle);
        float ax2 = a.x * a.x;

        return
        {
            c + ((1 - c) * a.x * a.x), ((1 - c) * a.x * a.y) + (s * a.z), ((1 - c) * a.x * a.z) - (s * a.y),
            ((1 - c) * a.x * a.y) - (s * a.z), c + ((1 - c) * a.y * a.y), ((1 - c) * a.y * a.z) + (s * a.x),
            ((1 - c) * a.x * a.z) + (s * a.y), ((1 - c) * a.y * a.z) - (s * a.x), c + ((1 - c) * a.z * a.z),
            0.0f, 0.0f, 0.0
        };
    }

    GEM_INLINE float4x3 GEM_VECTORCALL float4x3::rotate_axis_x(const float angle)
    {
        float s = sin(angle);
        float c = cos(angle);
        return
        {
            1.0f, 0.0f, 0.0f,
            0.0f, c, s,
            0.0f, -s, c,
            0.0f, 0.0f, 0.0
        };
    }

    GEM_INLINE  float4x3 GEM_VECTORCALL float4x3::rotate_axis_y(const float angle)
    {
        float s = sin(angle);
        float c = cos(angle);
        return
        {
            c, 0.0f, -s,
            0.0f, 1.0f, 0.0f,
            s, 0.0f, c,
            0.0f, 0.0f, 0.0
        };
    }

    GEM_INLINE float4x3 GEM_VECTORCALL float4x3::rotate_axis_z(const float angle)
    {
        float s = sin(angle);
        float c = cos(angle);
        return
        {
            c, s, 0.0f,
            -s, c, 0.0f,
            0.0f, 0.0f, 1.0f,
            0.0f, 0.0f, 0.0
        };
    }

    GEM_INLINE float4x3 GEM_VECTORCALL float4x3::rotate_forward_lh(const float3& direction, const float3& position, const float3& up)
    {
        float3 f = normalize(direction);
        float3 s = normalize(cross(up, f));
        float3 u = cross(f, s);

        return
        {
            s.x, s.y, s.z,
            u.x, u.y, u.z,
            f.x, f.y, f.z,
            position.x, position.y, position.z
        };
    }

    GEM_INLINE float4x3 GEM_VECTORCALL float4x3::rotate_forward_rh(const float3& direction, const float3& position, const float3& up)
    {
        return rotate_forward_lh(-direction, position, up);
    }

    GEM_INLINE float4x3 GEM_VECTORCALL float4x3::scale(const float3& s)
    {
        return
        {
            s.x, 0.0f, 0.0f,
            0.0f, s.y, 0.0f,
            0.0f, 0.0f, s.z,
            0.0f, 0.0f, 0.0
        };
    }

    GEM_INLINE float4x3 GEM_VECTORCALL float4x3::translate(const float3& t)
    {
        return
        {
            1.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 1.0f,
            t.x, t.y, t.z
        };
    }

    GEM_INLINE float4x3::float4x3(
        const float u1, const float u2, const float u3,
        const float v1, const float v2, const float v3,
        const float w1, const float w2, const float w3,
        const float t1, const float t2, const float t3)
        : m00(u1), m01(u2), m02(u3)
        , m10(v1), m11(v2), m12(v3)
        , m20(w1), m21(w2), m22(w3)
        , m30(t1), m31(t2), m32(t3)
    {

    }

    GEM_INLINE float4x3::float4x3(const float* o)
        : m00(o[0]), m01(o[1]), m02(o[2])
        , m10(o[3]), m11(o[4]), m12(o[5])
        , m20(o[6]), m21(o[7]), m22(o[8])
        , m30(o[9]), m31(o[10]), m32(o[11])
    {

    }

    GEM_INLINE float4x3::float4x3(const float3& u, const float3& v, const float3& w, const float3& t)
        : m00(u.x), m01(u.y), m02(u.z)
        , m10(v.x), m11(v.y), m12(v.z)
        , m20(w.x), m21(w.y), m22(w.z)
        , m30(t.x), m31(t.y), m32(t.z)
    {

    };

    GEM_INLINE float4x3::float4x3(const float4x3& o)
        : m00(o.m00), m01(o.m01), m02(o.m02)
        , m10(o.m10), m11(o.m11), m12(o.m12)
        , m20(o.m20), m21(o.m21), m22(o.m22)
        , m30(o.m30), m31(o.m31), m32(o.m32)
    {

    }

    GEM_INLINE float4x3 float4x3::transpose() const
    {
        return
        {
            m00, m10, m20, m30,
            m01, m11, m12, m31,
            m02, m12, m22, m32
        };
    }

    GEM_INLINE float4x3 float4x3::inverse() const
    {
        float4x3 out;
        float c11 = +((entries1d[4] * entries1d[8]) - (entries1d[5] * entries1d[7]));
        float c12 = -((entries1d[3] * entries1d[8]) - (entries1d[5] * entries1d[6]));
        float c13 = +((entries1d[3] * entries1d[7]) - (entries1d[4] * entries1d[6]));
        float c14 = -((((entries1d[4] * entries1d[8]) - (entries1d[5] * entries1d[7])) * entries1d[9]) + (((entries1d[5] * entries1d[6]) - (entries1d[3] * entries1d[8])) * entries1d[10]) + (((entries1d[3] * entries1d[7]) - (entries1d[4] * entries1d[6])) * entries1d[11]));

        float det = (entries1d[0] * c11) + (entries1d[1] * c12) + (entries1d[2] * c13);

        float c21 = -((entries1d[1] * entries1d[8]) - (entries1d[2] * entries1d[7]));
        float c22 = +((entries1d[0] * entries1d[8]) - (entries1d[2] * entries1d[6]));
        float c23 = -((entries1d[0] * entries1d[7]) - (entries1d[1] * entries1d[6]));
        float c24 = +((((entries1d[1] * entries1d[8]) - (entries1d[2] * entries1d[7])) * entries1d[9]) + (((entries1d[2] * entries1d[6]) - (entries1d[0] * entries1d[8])) * entries1d[10]) + (((entries1d[0] * entries1d[7]) - (entries1d[1] * entries1d[6])) * entries1d[11]));

        float c31 = +((entries1d[1] * entries1d[5]) - (entries1d[2] * entries1d[4]));
        float c32 = -((entries1d[0] * entries1d[5]) - (entries1d[2] * entries1d[3]));
        float c33 = +((entries1d[0] * entries1d[4]) - (entries1d[1] * entries1d[3]));
        float c34 = -((((entries1d[1] * entries1d[5]) - (entries1d[2] * entries1d[4])) * entries1d[9]) + (((entries1d[2] * entries1d[3]) - (entries1d[0] * entries1d[5])) * entries1d[10]) + (((entries1d[0] * entries1d[4]) - (entries1d[1] * entries1d[3])) * entries1d[11]));

        // float c41 = float(0);
        // float c42 = float(0);
        // float c43 = float(0);
        float c44 = +((((entries1d[1] * entries1d[5]) - (entries1d[2] * entries1d[4])) * entries1d[6]) + (((entries1d[2] * entries1d[3]) - (entries1d[0] * entries1d[5])) * entries1d[7]) + (((entries1d[0] * entries1d[4]) - (entries1d[1] * entries1d[3])) * entries1d[8]));

        float inv = 1.0f / det;

        out.entries1d[0] = c11 * inv;
        out.entries1d[1] = c21 * inv;
        out.entries1d[2] = c31 * inv;
        out.entries1d[3] = c12 * inv;
        out.entries1d[4] = c22 * inv;
        out.entries1d[5] = c32 * inv;
        out.entries1d[6] = c13 * inv;
        out.entries1d[7] = c23 * inv;
        out.entries1d[8] = c33 * inv;
        out.entries1d[9] = c14 * inv;
        out.entries1d[10] = c24 * inv;
        out.entries1d[11] = c34 * inv;
        return out;
    }

    GEM_INLINE float float4x3::determinant() const
    {
        float c11 = +((entries1d[4] * entries1d[8]) - (entries1d[5] * entries1d[7]));
        float c12 = -((entries1d[3] * entries1d[8]) - (entries1d[5] * entries1d[6]));
        float c13 = +((entries1d[3] * entries1d[7]) - (entries1d[4] * entries1d[6]));
        float c14 = -((((entries1d[4] * entries1d[8]) - (entries1d[5] * entries1d[7])) * entries1d[9]) + (((entries1d[5] * entries1d[6]) - (entries1d[3] * entries1d[8])) * entries1d[10]) + (((entries1d[3] * entries1d[7]) - (entries1d[4] * entries1d[6])) * entries1d[11]));

        return (entries1d[0] * c11) + (entries1d[1] * c12) + (entries1d[2] * c13);
    }

    GEM_INLINE float3 float4x3::euler() const
    {
        float3 out;
        float sp = -entries1d[7];
        if (sp <= 1.0f)
        {
            out.x = -1.570796f;
        }
        else if (sp >= 1.0f)
        {
            out.x = -1.570796f;
        }
        else
        {
            out.x = asin(sp);
        }

        if (abs(sp) > 0.9999f)
        {
            out.z = 0.0f;
            out.y = atan2(-entries1d[2], entries1d[0]);
        }
        else
        {
            out.y = atan2(entries1d[6], entries1d[8]);
            out.z = atan2(entries1d[1], entries1d[4]);
        }
        return out;
    }

    GEM_INLINE float4x3& GEM_VECTORCALL float4x3::operator=(const float4x3& rhs)
    {
        m00 = rhs.m00; m01 = rhs.m01; m02 = rhs.m02;
        m10 = rhs.m10; m11 = rhs.m11; m12 = rhs.m12;
        m20 = rhs.m20; m21 = rhs.m21; m22 = rhs.m22;
        m30 = rhs.m30; m31 = rhs.m31; m32 = rhs.m32;
        return *this;
    }

    GEM_INLINE float4x3 GEM_VECTORCALL float4x3::operator+=(const float4x3& rhs)
    {
        m00 += rhs.m00; m01 += rhs.m01; m02 += rhs.m02;
        m10 += rhs.m10; m11 += rhs.m11; m12 += rhs.m12;
        m20 += rhs.m20; m21 += rhs.m21; m22 += rhs.m22;
        m30 += rhs.m30; m31 += rhs.m31; m32 += rhs.m32;
        return *this;
    }

    GEM_INLINE float4x3 GEM_VECTORCALL float4x3::operator-=(const float4x3& rhs)
    {
        m00 -= rhs.m00; m01 -= rhs.m01; m02 -= rhs.m02;
        m10 -= rhs.m10; m11 -= rhs.m11; m12 -= rhs.m12;
        m20 -= rhs.m20; m21 -= rhs.m21; m22 -= rhs.m22;
        m30 -= rhs.m30; m31 -= rhs.m31; m32 -= rhs.m32;
        return *this;
    }

    GEM_INLINE float4x3 GEM_VECTORCALL float4x3::operator*=(const float rhs)
    {
        m00 *= rhs; m01 *= rhs; m02 *= rhs;
        m10 *= rhs; m11 *= rhs; m12 *= rhs;
        m20 *= rhs; m21 *= rhs; m22 *= rhs;
        m30 *= rhs; m31 *= rhs; m32 *= rhs;
        return *this;
    }

    GEM_INLINE float3& GEM_VECTORCALL float4x3::operator[](const unsigned int index)
    {
        return rows[index];
    }

    GEM_INLINE float3 GEM_VECTORCALL float4x3::operator[](const unsigned int index) const
    {
        return rows[index];
    }

    GEM_INLINE float GEM_VECTORCALL float4x3::operator()(const unsigned int row, const unsigned int column)
    {
        return entries2d[row][column];
    }

    GEM_INLINE float4x3 float4x3::operator-() const
    {
        return
        {
            -m00, -m01, -m02,
            -m10, -m11, -m12,
            -m20, -m21, -m22,
            -m30, -m31, -m32
        };
    }

    GEM_INLINE float4x3 GEM_VECTORCALL operator+(const float4x3& lhs, const float4x3& rhs)
    {
        return
        {
            lhs.entries1d[0] + rhs.entries1d[0], lhs.entries1d[1] + rhs.entries1d[1], lhs.entries1d[2] + rhs.entries1d[2],
            lhs.entries1d[3] + rhs.entries1d[3], lhs.entries1d[4] + rhs.entries1d[4], lhs.entries1d[5] + rhs.entries1d[5],
            lhs.entries1d[6] + rhs.entries1d[6], lhs.entries1d[7] + rhs.entries1d[7], lhs.entries1d[8] + rhs.entries1d[8],
            lhs.entries1d[9] + rhs.entries1d[9], lhs.entries1d[10] + rhs.entries1d[10], lhs.entries1d[11] + rhs.entries1d[11]
        };
    }

    GEM_INLINE float4x3 GEM_VECTORCALL operator-(const float4x3& lhs, const float4x3& rhs)
    {
        return
        {
            lhs.entries1d[0] - rhs.entries1d[0], lhs.entries1d[1] - rhs.entries1d[1], lhs.entries1d[2] - rhs.entries1d[2],
            lhs.entries1d[3] - rhs.entries1d[3], lhs.entries1d[4] - rhs.entries1d[4], lhs.entries1d[5] - rhs.entries1d[5],
            lhs.entries1d[6] - rhs.entries1d[6], lhs.entries1d[7] - rhs.entries1d[7], lhs.entries1d[8] - rhs.entries1d[8],
            lhs.entries1d[9] - rhs.entries1d[9], lhs.entries1d[10] - rhs.entries1d[10], lhs.entries1d[11] - rhs.entries1d[11]
        };
    }

    GEM_INLINE float4x3 GEM_VECTORCALL operator*(const float4x3& lhs, const float4x3& rhs)
    {
        return
        {
            (lhs.entries1d[0] * rhs.entries1d[0]) + (lhs.entries1d[1] * rhs.entries1d[3]) + (lhs.entries1d[2] * rhs.entries1d[6]),
            (lhs.entries1d[0] * rhs.entries1d[1]) + (lhs.entries1d[1] * rhs.entries1d[4]) + (lhs.entries1d[2] * rhs.entries1d[7]),
            (lhs.entries1d[0] * rhs.entries1d[2]) + (lhs.entries1d[1] * rhs.entries1d[5]) + (lhs.entries1d[2] * rhs.entries1d[8]),
            (lhs.entries1d[3] * rhs.entries1d[0]) + (lhs.entries1d[4] * rhs.entries1d[3]) + (lhs.entries1d[5] * rhs.entries1d[6]),
            (lhs.entries1d[3] * rhs.entries1d[1]) + (lhs.entries1d[4] * rhs.entries1d[4]) + (lhs.entries1d[5] * rhs.entries1d[7]),
            (lhs.entries1d[3] * rhs.entries1d[2]) + (lhs.entries1d[4] * rhs.entries1d[5]) + (lhs.entries1d[5] * rhs.entries1d[8]),
            (lhs.entries1d[6] * rhs.entries1d[0]) + (lhs.entries1d[7] * rhs.entries1d[3]) + (lhs.entries1d[8] * rhs.entries1d[6]),
            (lhs.entries1d[6] * rhs.entries1d[1]) + (lhs.entries1d[7] * rhs.entries1d[4]) + (lhs.entries1d[8] * rhs.entries1d[7]),
            (lhs.entries1d[6] * rhs.entries1d[2]) + (lhs.entries1d[7] * rhs.entries1d[5]) + (lhs.entries1d[8] * rhs.entries1d[8]),
            (lhs.entries1d[9] * rhs.entries1d[0]) + (lhs.entries1d[10] * rhs.entries1d[3]) + (lhs.entries1d[11] * rhs.entries1d[6]) + rhs.entries1d[9],
            (lhs.entries1d[9] * rhs.entries1d[1]) + (lhs.entries1d[10] * rhs.entries1d[4]) + (lhs.entries1d[11] * rhs.entries1d[7]) + rhs.entries1d[10],
            (lhs.entries1d[9] * rhs.entries1d[2]) + (lhs.entries1d[10] * rhs.entries1d[5]) + (lhs.entries1d[11] * rhs.entries1d[8]) + rhs.entries1d[11]
        };
    }

    GEM_INLINE float4x3 GEM_VECTORCALL operator*(const float lhs, const float4x3& rhs)
    {
        return
        {
            rhs.entries1d[0] * lhs, rhs.entries1d[ 1] * lhs, rhs.entries1d[ 2] * lhs,
            rhs.entries1d[3] * lhs, rhs.entries1d[ 4] * lhs, rhs.entries1d[ 5] * lhs,
            rhs.entries1d[6] * lhs, rhs.entries1d[ 7] * lhs, rhs.entries1d[ 8] * lhs,
            rhs.entries1d[9] * lhs, rhs.entries1d[10] * lhs, rhs.entries1d[11] * lhs
        };
    }

    GEM_INLINE float4x3 GEM_VECTORCALL operator*(const float4x3& lhs, const float rhs)
    {
        return
        {
            lhs.entries1d[0] * rhs, lhs.entries1d[1] * rhs, lhs.entries1d[2] * rhs,
            lhs.entries1d[3] * rhs, lhs.entries1d[4] * rhs, lhs.entries1d[5] * rhs,
            lhs.entries1d[6] * rhs, lhs.entries1d[7] * rhs, lhs.entries1d[8] * rhs,
            lhs.entries1d[9] * rhs, lhs.entries1d[10] * rhs, lhs.entries1d[11] * rhs
        };
    }

    GEM_INLINE float3 GEM_VECTORCALL operator*(const float3& lhs, const float4x3& rhs)
    {
        return
        {
            (lhs.components[0] * rhs.entries1d[0]) + (lhs.components[1] * rhs.entries1d[3]) + (lhs.components[2] * rhs.entries1d[6]) + rhs.entries1d[9],
            (lhs.components[0] * rhs.entries1d[1]) + (lhs.components[1] * rhs.entries1d[4]) + (lhs.components[2] * rhs.entries1d[7]) + rhs.entries1d[10],
            (lhs.components[0] * rhs.entries1d[2]) + (lhs.components[1] * rhs.entries1d[5]) + (lhs.components[2] * rhs.entries1d[8]) + rhs.entries1d[11]
        };
    }

#pragma endregion

#pragma region float4x4

    struct float4x4
    {
        union
        {
            struct
            {
                float m00, m01, m02, m03;
                float m10, m11, m12, m13;
                float m20, m21, m22, m23;
                float m30, m31, m32, m33;
            };

            struct
            {
                float4 u, v, w, t;
            };

            float4 rows[4];
            float entries1d[16];
            float entries2d[4][4];
        };

        static float4x4 identity();

        static float4x4 GEM_VECTORCALL rotate_euler(const float pitch, const float yaw, const float roll);

        static float4x4 GEM_VECTORCALL rotate_axis_angle(const float3& a, const float angle);

        static float4x4 GEM_VECTORCALL rotate_axis_x(const float angle);

        static float4x4 GEM_VECTORCALL rotate_axis_y(const float angle);

        static float4x4 GEM_VECTORCALL rotate_axis_z(const float angle);

        static float4x4 GEM_VECTORCALL rotate_forward_lh(const float3& direction, const float3& position, const float3& up);

        static float4x4 GEM_VECTORCALL rotate_forward_rh(const float3& direction, const float3& position, const float3& up);

        static float4x4 GEM_VECTORCALL scale(const float3& s);

        static float4x4 GEM_VECTORCALL translate(const float3& t);

        static float4x4 GEM_VECTORCALL orthographic_lh(const float left, const float right, const float bottom, const float top, const float zNear, const float zFar);

        static float4x4 GEM_VECTORCALL orthographic_rh(const float left, const float right, const float bottom, const float top, const float zNear, const float zFar);

        static float4x4 GEM_VECTORCALL orthographic_normalized_lh(const float left, const float right, const float bottom, const float top, const float zNear, const float zFar);

        static float4x4 GEM_VECTORCALL orthographic_normalized_rh(const float left, const float right, const float bottom, const float top, const float zNear, const float zFar);

        static float4x4 GEM_VECTORCALL perspective_lh(const float fovy, const float aspectRatio, const float zNear, const float zFar);

        static float4x4 GEM_VECTORCALL perspective_rh(const float fovy, const float aspectRatio, const float zNear, const float zFar);

        static float4x4 GEM_VECTORCALL perspective_normalized_lh(const float fovy, const float aspectRatio, const float zNear, const float zFar);

        static float4x4 GEM_VECTORCALL perspective_normalized_rh(const float fovy, const float aspectRatio, const float zNear, const float zFar);

        static float4x4 GEM_VECTORCALL look_at_lh(const float3& target, const float3& position, const float3& up);

        static float4x4 GEM_VECTORCALL look_at_rh(const float3& target, const float3& position, const float3& up);

        static float4x4 GEM_VECTORCALL look_to_lh(const float3& direction, const float3& position, const float3& up);

        static float4x4 GEM_VECTORCALL look_to_rh(const float3& direction, const float3& position, const float3& up);

        float4x4(const float u1, const float u2, const float u3, const float u4,
            const float v1, const float v2, const float v3, const float v4,
            const float w1, const float w2, const float w3, const float w4,
            const float t1, const float t2, const float t3, const float t4);

        float4x4(const float4& u, const float4& v, const float4& w, const float4& t);

        float4x4(const float4x4& other);

        float4x4(const float* o);

        float4x4() = default;

        float4x4 transpose() const;

        float4x4 inverse() const;

        float determinant()	const;

        float3 euler() const;

        float4x4& GEM_VECTORCALL operator=(const float4x4& rhs);

        float4x4 GEM_VECTORCALL operator+=(const float4x4& rhs);

        float4x4 GEM_VECTORCALL operator-=(const float4x4& rhs);

        float4x4 GEM_VECTORCALL operator*=(const float rhs);

        float4x4 GEM_VECTORCALL operator-() const;

        float4& GEM_VECTORCALL operator[](const unsigned int index);

        float4 GEM_VECTORCALL operator[](const unsigned int index) const;

        float GEM_VECTORCALL operator()(const unsigned int row, const unsigned int column);
    };

    float4x4 GEM_VECTORCALL operator+(const float4x4& lhs, const float4x4& rhs);

    float4x4 GEM_VECTORCALL operator-(const float4x4& lhs, const float4x4& rhs);

    float4x4 GEM_VECTORCALL operator*(const float4x4& lhs, const float4x4& rhs);

    float4x4 GEM_VECTORCALL operator*(const float4x4& lhs, const float rhs);

    float4x4 GEM_VECTORCALL operator*(const float lhs, const float4x4& rhs);

    float4 GEM_VECTORCALL operator*(const float4& lhs, const float4x4& rhs);

    GEM_INLINE float4x4 float4x4::identity()
    {
        return
        {
            1.0f, 0.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 0.0f, 1.0f
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::rotate_euler(const float pitch, const float yaw, const float roll)
    {
        float ch = cosf(yaw); float cp = cosf(pitch); float cb = cosf(roll);
        float sh = sinf(yaw); float sp = sinf(pitch); float sb = sinf(roll);

        return
        {
             ch * cb + sh * sp * sb, sb * cp, -sh * cb + ch * sp * sb, 0.0f,
            -ch * sb + sh * sp * cb, cb * cp,  sb * sh + ch * sp * cb, 0.0f,
             sh * cp, -sp,  ch * cp, 0.0f,
             0.0f, 0.0f, 0.0f, 1.0f
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::rotate_axis_angle(const float3& a, const float angle)
    {
        float s = sin(angle);
        float c = cos(angle);

        return
        {
            c + ((1 - c) * a.x * a.x), ((1 - c) * a.x * a.y) + (s * a.z), ((1 - c) * a.x * a.z) - (s * a.y), 0.0f,
            ((1 - c) * a.x * a.y) - (s * a.z), c + ((1 - c) * a.y * a.y), ((1 - c) * a.y * a.z) + (s * a.x), 0.0f,
            ((1 - c) * a.x * a.z) + (s * a.y), ((1 - c) * a.y * a.z) - (s * a.x), c + ((1 - c) * a.z * a.z), 0.0f,
            0.0f, 0.0f, 0.0, 1.0f
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::rotate_axis_x(const float angle)
    {
        float s = sin(angle);
        float c = cos(angle);
        return
        {
            1.0f, 0.0f, 0.0f, 0.0f,
            0.0f, c, s, 0.0f,
            0.0f, -s, c, 0.0f,
            0.0f, 0.0f, 0.0, 1.0f
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::rotate_axis_y(const float angle)
    {
        float s = sin(angle);
        float c = cos(angle);
        return
        {
            c, 0.0f, -s, 0.0f,
            0.0f, 1.0f, 0.0f, 0.0f,
            s, 0.0f, c, 0.0f,
            0.0f, 0.0f, 0.0, 1.0f
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::rotate_axis_z(const float angle)
    {
        float s = sin(angle);
        float c = cos(angle);
        return
        {
            c, s, 0.0f, 0.0f,
            -s, c, 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 0.0, 1.0f
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::rotate_forward_lh(const float3& direction, const float3& position, const float3& up)
    {
        float3 f = normalize(direction);
        float3 s = normalize(cross(up, f));
        float3 u = cross(f, s);

        return
        {
            s.x, s.y, s.z, 0.0f,
            u.x, u.y, u.z, 0.0f,
            f.x, f.y, f.z, 0.0f,
            position.x, position.y, position.z, 1.0f
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::rotate_forward_rh(const float3& direction, const float3& position, const float3& up)
    {
        return rotate_forward_lh(-direction, position, up);
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::scale(const float3& s)
    {
        return
        {
            s.x, 0.0f, 0.0f, 0.0f,
            0.0f, s.y, 0.0f, 0.0f,
            0.0f, 0.0f, s.z, 0.0f,
            0.0f, 0.0f, 0.0, 1.0f
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::translate(const float3& t)
    {
        return
        {
            1.0f, 0.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f, 0.0f,
            t.x, t.y, t.z, 1.0f
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::orthographic_lh(
        const float left, const float right,
        const float bottom, const float top,
        const float zNear, const float zFar)
    {
        return
        {
            2.0f / (right - left), 0.0f, 0.0f, 0.0f,
            0.0f, 2.0f / (top - bottom), 0.0f, 0.0f,
            0.0f, 0.0f, 2.0f / (zFar - zNear), 0.0f,
            0.0f, 0.0f, -((zFar + zNear) / (zFar - zNear)), 1.0f
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::orthographic_rh(
        const float left, const float right,
        const float bottom, const float top,
        const float zNear, const float zFar)
    {
        return
        {
            2.0f / (right - left), 0.0f, 0.0f, 0.0f,
            0.0f, 2.0f / (top - bottom), 0.0f, 0.0f,
            0.0f, 0.0f, -2.0f / (zFar - zNear), 0.0f,
            0.0f, 0.0f, -((zFar + zNear) / (zFar - zNear)), 1.0f
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::orthographic_normalized_lh(
        const float left, const float right,
        const float bottom, const float top,
        const float zNear, const float zFar)
    {
        return
        {
            2.0f / (right - left), 0.0f, 0.0f, 0.0f,
            0.0f, 2.0f / (top - bottom), 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f / (zFar - zNear), 0.0f,
            0.0f, 0.0f, zNear / (zNear - zFar), 1.0f
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::orthographic_normalized_rh(
        const float left, const float right,
        const float bottom, const float top,
        const float zNear, const float zFar)
    {
        return
        {
            2.0f / (right - left), 0.0f, 0.0f, 0.0f,
            0.0f, 2.0f / (top - bottom), 0.0f, 0.0f,
            0.0f, 0.0f, -1.0f / (zFar - zNear), 0.0f,
            0.0f, 0.0f, -(zNear / (zNear - zFar)), 1.0f
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::perspective_lh(
        const float fovy, const float aspectRatio,
        const float zNear, const float zFar)
    {
        float tanHalfFovy = tanf(fovy * 0.5f);
        float zoomY = 1.0f / tanHalfFovy;
        float zoomX = 1.0f / (aspectRatio * tanHalfFovy);

        return
        {
            zoomX, 0.0f, 0.0f, 0.0f,
            0.0f, zoomY, 0.0f, 0.0f,
            0.0f, 0.0f, ((zFar + zNear) / (zFar - zNear)), 1.0f,
            0.0f, 0.0f, (-2.0f * zNear * zFar) / (zFar - zNear), 0.0f
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::perspective_rh(
        const float fovy, const float aspectRatio,
        const float zNear, const float zFar)
    {
        float tanHalfFovy = tanf(fovy * 0.5f);
        float zoomY = 1.0f / tanHalfFovy;
        float zoomX = 1.0f / (aspectRatio * tanHalfFovy);

        return
        {
            zoomX, 0.0f, 0.0f, 0.0f,
            0.0f, zoomY, 0.0f, 0.0f,
            0.0f, 0.0f, -((zFar + zNear) / (zFar - zNear)), -1.0f,
            0.0f, 0.0f, (-2.0f * zNear * zFar) / (zFar - zNear), 0.0f
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::perspective_normalized_lh(
        const float fovy, const float aspectRatio, const float zNear, const float zFar)
    {
        float tanHalfFovy = tanf(fovy * 0.5f);
        float zoomY = 1.0f / tanHalfFovy;
        float zoomX = 1.0f / (aspectRatio * tanHalfFovy);

        return
        {
            zoomX, 0.0f, 0.0f, 0.0f,
            0.0f, zoomY, 0.0f, 0.0f,
            0.0f, 0.0f, zFar / (zFar - zNear), 1.0f,
            0.0f, 0.0f, (-zNear * zFar) / (zFar - zNear), 0.0f
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::perspective_normalized_rh(
        const float fovy, const float aspectRatio, const float zNear, const float zFar)
    {
        float tanHalfFovy = tanf(fovy * 0.5f);
        float zoomY = 1.0f / tanHalfFovy;
        float zoomX = 1.0f / (aspectRatio * tanHalfFovy);

        return
        {
            zoomX, 0.0f, 0.0f, 0.0f,
            0.0f, zoomY, 0.0f, 0.0f,
            0.0f, 0.0f, -zFar / (zFar - zNear), -1.0f,
            0.0f, 0.0f, (-zNear * zFar) / (zFar - zNear), 0.0f
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::look_to_lh(const float3& direction, const float3& position, const float3& up)
    {
        return float4x4::look_at_lh(position + direction, position, up);
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::look_to_rh(const float3& direction, const float3& position, const float3& up)
    {
        return float4x4::look_at_rh(position + direction, position, up);
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::look_at_lh(const float3& target, const float3& position, const float3& up)
    {
        float3 f = normalize(target - position);
        float3 s = normalize(cross(up, f));
        float3 u = cross(f, s);

        return
        {
            s.x, u.x, f.x, 0.0f,
            s.y, u.y, f.y, 0.0f,
            s.z, u.z, f.z, 0.0f,
            -dot(s, position), -dot(u, position), -dot(f, position), 1.0f,
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::look_at_rh(const float3& target, const float3& position, const float3& up)
    {
        float3 f = normalize(position - target);
        float3 s = normalize(cross(up, f));
        float3 u = cross(f, s);

        return
        {
            s.x, u.x, f.x, 0.0f,
            s.y, u.y, f.y, 0.0f,
            s.z, u.z, f.z, 0.0f,
            -dot(s, position), -dot(u, position), -dot(f, position), 1.0f,
        };
    }

    GEM_INLINE float4x4::float4x4(
        const float u1, const float u2, const float u3, const float u4,
        const float v1, const float v2, const float v3, const float v4,
        const float w1, const float w2, const float w3, const float w4,
        const float t1, const float t2, const float t3, const float t4)
        : m00(u1), m01(u2), m02(u3), m03(u4)
        , m10(v1), m11(v2), m12(v3), m13(v4)
        , m20(w1), m21(w2), m22(w3), m23(w4)
        , m30(t1), m31(t2), m32(t3), m33(t4)

    {

    }

    GEM_INLINE float4x4::float4x4(const float* o)
        : m00(o[ 0]), m01(o[ 1]), m02(o[ 2]), m03(o[ 3])
        , m10(o[ 4]), m11(o[ 5]), m12(o[ 6]), m13(o[ 7])
        , m20(o[ 8]), m21(o[ 9]), m22(o[10]), m23(o[11])
        , m30(o[12]), m31(o[13]), m32(o[14]), m33(o[15])
    {

    }

    GEM_INLINE float4x4::float4x4(const float4& u, const float4& v, const float4& w, const float4& t)
        : m00(u.x), m01(u.y), m02(u.z), m03(u.w)
        , m10(v.x), m11(v.y), m12(v.z), m13(v.w)
        , m20(w.x), m21(w.y), m22(w.z), m23(w.w)
        , m30(t.x), m31(t.y), m32(t.z), m33(t.w)
    {

    };

    GEM_INLINE float4x4::float4x4(const float4x4& o)
        : m00(o.m00), m01(o.m01), m02(o.m02), m03(o.m03)
        , m10(o.m10), m11(o.m11), m12(o.m12), m13(o.m13)
        , m20(o.m20), m21(o.m21), m22(o.m22), m23(o.m23)
        , m30(o.m30), m31(o.m31), m32(o.m32), m33(o.m33)
    {

    }

    GEM_INLINE float4x4 float4x4::transpose() const
    {
        return
        {
            m00, m10, m20, m30,
            m01, m11, m21, m31,
            m02, m12, m22, m32,
            m03, m13, m23, m33
        };
    }

    GEM_INLINE float4x4 float4x4::inverse() const
    {
        float4x4 out;

        float c11 = +((((entries1d[6] * entries1d[11]) - (entries1d[7] * entries1d[10])) * entries1d[13]) + (((entries1d[7] * entries1d[9]) - (entries1d[5] * entries1d[11])) * entries1d[14]) + (((entries1d[5] * entries1d[10]) - (entries1d[6] * entries1d[9])) * entries1d[15]));
        float c12 = -((((entries1d[6] * entries1d[11]) - (entries1d[7] * entries1d[10])) * entries1d[12]) + (((entries1d[7] * entries1d[8]) - (entries1d[4] * entries1d[11])) * entries1d[14]) + (((entries1d[4] * entries1d[10]) - (entries1d[6] * entries1d[8])) * entries1d[15]));
        float c13 = +((((entries1d[5] * entries1d[11]) - (entries1d[7] * entries1d[9])) * entries1d[12]) + (((entries1d[7] * entries1d[8]) - (entries1d[4] * entries1d[11])) * entries1d[13]) + (((entries1d[4] * entries1d[9]) - (entries1d[5] * entries1d[8])) * entries1d[15]));
        float c14 = -((((entries1d[5] * entries1d[10]) - (entries1d[6] * entries1d[9])) * entries1d[12]) + (((entries1d[6] * entries1d[8]) - (entries1d[4] * entries1d[10])) * entries1d[13]) + (((entries1d[4] * entries1d[9]) - (entries1d[5] * entries1d[8])) * entries1d[14]));

        float det = (entries1d[0] * c11) + (entries1d[1] * c12) + (entries1d[2] * c13) + (entries1d[3] * c14);

        float c21 = -((((entries1d[2] * entries1d[11]) - (entries1d[3] * entries1d[10])) * entries1d[13]) + (((entries1d[3] * entries1d[9]) - (entries1d[1] * entries1d[11])) * entries1d[14]) + (((entries1d[1] * entries1d[10]) - (entries1d[2] * entries1d[9])) * entries1d[15]));
        float c22 = +((((entries1d[2] * entries1d[11]) - (entries1d[3] * entries1d[10])) * entries1d[12]) + (((entries1d[3] * entries1d[8]) - (entries1d[0] * entries1d[11])) * entries1d[14]) + (((entries1d[0] * entries1d[10]) - (entries1d[2] * entries1d[8])) * entries1d[15]));
        float c23 = -((((entries1d[1] * entries1d[11]) - (entries1d[3] * entries1d[9])) * entries1d[12]) + (((entries1d[3] * entries1d[8]) - (entries1d[0] * entries1d[11])) * entries1d[13]) + (((entries1d[0] * entries1d[9]) - (entries1d[1] * entries1d[8])) * entries1d[15]));
        float c24 = +((((entries1d[1] * entries1d[10]) - (entries1d[2] * entries1d[9])) * entries1d[12]) + (((entries1d[2] * entries1d[8]) - (entries1d[0] * entries1d[10])) * entries1d[13]) + (((entries1d[0] * entries1d[9]) - (entries1d[1] * entries1d[8])) * entries1d[14]));

        float c31 = +((((entries1d[2] * entries1d[7]) - (entries1d[3] * entries1d[6])) * entries1d[13]) + (((entries1d[3] * entries1d[5]) - (entries1d[1] * entries1d[7])) * entries1d[14]) + (((entries1d[1] * entries1d[6]) - (entries1d[2] * entries1d[5])) * entries1d[15]));
        float c32 = -((((entries1d[2] * entries1d[7]) - (entries1d[3] * entries1d[6])) * entries1d[12]) + (((entries1d[3] * entries1d[4]) - (entries1d[0] * entries1d[7])) * entries1d[14]) + (((entries1d[0] * entries1d[6]) - (entries1d[2] * entries1d[4])) * entries1d[15]));
        float c33 = +((((entries1d[1] * entries1d[7]) - (entries1d[3] * entries1d[5])) * entries1d[12]) + (((entries1d[3] * entries1d[4]) - (entries1d[0] * entries1d[7])) * entries1d[13]) + (((entries1d[0] * entries1d[5]) - (entries1d[1] * entries1d[4])) * entries1d[15]));
        float c34 = -((((entries1d[1] * entries1d[6]) - (entries1d[2] * entries1d[5])) * entries1d[12]) + (((entries1d[2] * entries1d[4]) - (entries1d[0] * entries1d[6])) * entries1d[13]) + (((entries1d[0] * entries1d[5]) - (entries1d[1] * entries1d[4])) * entries1d[14]));

        float c41 = -((((entries1d[2] * entries1d[7]) - (entries1d[3] * entries1d[6])) * entries1d[9]) + (((entries1d[3] * entries1d[5]) - (entries1d[1] * entries1d[7])) * entries1d[10]) + (((entries1d[1] * entries1d[6]) - (entries1d[2] * entries1d[5])) * entries1d[11]));
        float c42 = +((((entries1d[2] * entries1d[7]) - (entries1d[3] * entries1d[6])) * entries1d[8]) + (((entries1d[3] * entries1d[4]) - (entries1d[0] * entries1d[7])) * entries1d[10]) + (((entries1d[0] * entries1d[6]) - (entries1d[2] * entries1d[4])) * entries1d[11]));
        float c43 = -((((entries1d[1] * entries1d[7]) - (entries1d[3] * entries1d[5])) * entries1d[8]) + (((entries1d[3] * entries1d[4]) - (entries1d[0] * entries1d[7])) * entries1d[9]) + (((entries1d[0] * entries1d[5]) - (entries1d[1] * entries1d[4])) * entries1d[11]));
        float c44 = +((((entries1d[1] * entries1d[6]) - (entries1d[2] * entries1d[5])) * entries1d[8]) + (((entries1d[2] * entries1d[4]) - (entries1d[0] * entries1d[6])) * entries1d[9]) + (((entries1d[0] * entries1d[5]) - (entries1d[1] * entries1d[4])) * entries1d[10]));

        float inv = 1.0f / det;

        out.entries1d[0] = c11 * inv;
        out.entries1d[1] = c21 * inv;
        out.entries1d[2] = c31 * inv;
        out.entries1d[3] = c41 * inv;
        out.entries1d[4] = c12 * inv;
        out.entries1d[5] = c22 * inv;
        out.entries1d[6] = c32 * inv;
        out.entries1d[7] = c42 * inv;
        out.entries1d[8] = c13 * inv;
        out.entries1d[9] = c23 * inv;
        out.entries1d[10] = c33 * inv;
        out.entries1d[11] = c43 * inv;
        out.entries1d[12] = c14 * inv;
        out.entries1d[13] = c24 * inv;
        out.entries1d[14] = c34 * inv;
        out.entries1d[15] = c44 * inv;

        return out;
    }

    GEM_INLINE float float4x4::determinant() const
    {
        float c11 = +((((entries1d[6] * entries1d[11]) - (entries1d[7] * entries1d[10])) * entries1d[13]) + (((entries1d[7] * entries1d[9]) - (entries1d[5] * entries1d[11])) * entries1d[14]) + (((entries1d[5] * entries1d[10]) - (entries1d[6] * entries1d[9])) * entries1d[15]));
        float c12 = -((((entries1d[6] * entries1d[11]) - (entries1d[7] * entries1d[10])) * entries1d[12]) + (((entries1d[7] * entries1d[8]) - (entries1d[4] * entries1d[11])) * entries1d[14]) + (((entries1d[4] * entries1d[10]) - (entries1d[6] * entries1d[8])) * entries1d[15]));
        float c13 = +((((entries1d[5] * entries1d[11]) - (entries1d[7] * entries1d[9])) * entries1d[12]) + (((entries1d[7] * entries1d[8]) - (entries1d[4] * entries1d[11])) * entries1d[13]) + (((entries1d[4] * entries1d[9]) - (entries1d[5] * entries1d[8])) * entries1d[15]));
        float c14 = -((((entries1d[5] * entries1d[10]) - (entries1d[6] * entries1d[9])) * entries1d[12]) + (((entries1d[6] * entries1d[8]) - (entries1d[4] * entries1d[10])) * entries1d[13]) + (((entries1d[4] * entries1d[9]) - (entries1d[5] * entries1d[8])) * entries1d[14]));

        return (entries1d[0] * c11) + (entries1d[1] * c12) + (entries1d[2] * c13) + (entries1d[3] * c14);
    }

    GEM_INLINE float3 float4x4::euler() const
    {
        float3 out;
        float sp = -entries1d[9];
        if (sp <= 1.0f)
        {
            out.x = -1.690998f;
        }
        else if (sp >= 1.0f)
        {
            out.x = -1.690998f;
        }
        else
        {
            out.x = asin(sp);
        }

        if (abs(sp) > 0.9999f)
        {
            out.z = 0.0f;
            out.y = atan2(-entries1d[2], entries1d[0]);
        }
        else
        {
            out.y = atan2(entries1d[8], entries1d[10]);
            out.z = atan2(entries1d[1], entries1d[5]);
        }
        return out;
    }

    GEM_INLINE float4x4& GEM_VECTORCALL float4x4::operator=(const float4x4& rhs)
    {
        m00 = rhs.m00; m01 = rhs.m01; m02 = rhs.m02; m03 = rhs.m03;
        m10 = rhs.m10; m11 = rhs.m11; m12 = rhs.m12; m13 = rhs.m13;
        m20 = rhs.m20; m21 = rhs.m21; m22 = rhs.m22; m23 = rhs.m23;
        m30 = rhs.m30; m31 = rhs.m31; m32 = rhs.m32; m33 = rhs.m33;
        return *this;
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::operator+=(const float4x4& rhs)
    {
        m00 += rhs.m00; m01 += rhs.m01; m02 += rhs.m02; m03 += rhs.m03;
        m10 += rhs.m10; m11 += rhs.m11; m12 += rhs.m12; m13 += rhs.m13;
        m20 += rhs.m20; m21 += rhs.m21; m22 += rhs.m22; m23 += rhs.m23;
        m30 += rhs.m30; m31 += rhs.m31; m32 += rhs.m32; m33 += rhs.m33;
        return *this;
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::operator-=(const float4x4& rhs)
    {
        m00 -= rhs.m00; m01 -= rhs.m01; m02 -= rhs.m02; m03 -= rhs.m03;
        m10 -= rhs.m10; m11 -= rhs.m11; m12 -= rhs.m12; m13 -= rhs.m13;
        m20 -= rhs.m20; m21 -= rhs.m21; m22 -= rhs.m22; m23 -= rhs.m23;
        m30 -= rhs.m30; m31 -= rhs.m31; m32 -= rhs.m32; m33 -= rhs.m33;
        return *this;
    }

    GEM_INLINE float4x4 GEM_VECTORCALL float4x4::operator*=(const float rhs)
    {
        m00 *= rhs; m01 *= rhs; m02 *= rhs; m03 *= rhs;
        m10 *= rhs; m11 *= rhs; m12 *= rhs; m13 *= rhs;
        m20 *= rhs; m21 *= rhs; m22 *= rhs; m23 *= rhs;
        m30 *= rhs; m31 *= rhs; m32 *= rhs; m33 *= rhs;
        return *this;
    }

    GEM_INLINE float4x4 float4x4::operator-() const
    {
        return
        {
            -m00, -m01, -m02, -m03,
            -m10, -m11, -m12, -m13,
            -m20, -m21, -m22, -m23,
            -m30, -m31, -m32, -m33
        };
    }

    GEM_INLINE float4& GEM_VECTORCALL float4x4::operator[](const unsigned int index)
    {
        return rows[index];
    }

    GEM_INLINE float4 GEM_VECTORCALL float4x4::operator[](const unsigned int index) const
    {
        return rows[index];
    }

    GEM_INLINE float GEM_VECTORCALL float4x4::operator()(const unsigned int row, const unsigned int column)
    {
        return entries2d[row][column];
    }

    GEM_INLINE float4x4 GEM_VECTORCALL operator+(const float4x4& lhs, const float4x4& rhs)
    {
        return
        {
            lhs.entries1d[0] + rhs.entries1d[0], lhs.entries1d[1] + rhs.entries1d[1], lhs.entries1d[2] + rhs.entries1d[2], lhs.entries1d[3] + rhs.entries1d[3],
            lhs.entries1d[4] + rhs.entries1d[4], lhs.entries1d[5] + rhs.entries1d[5], lhs.entries1d[6] + rhs.entries1d[6], lhs.entries1d[7] + rhs.entries1d[7],
            lhs.entries1d[8] + rhs.entries1d[8], lhs.entries1d[9] + rhs.entries1d[9], lhs.entries1d[10] + rhs.entries1d[10], lhs.entries1d[11] + rhs.entries1d[11],
            lhs.entries1d[12] + rhs.entries1d[12], lhs.entries1d[13] + rhs.entries1d[13], lhs.entries1d[14] + rhs.entries1d[14], lhs.entries1d[15] + rhs.entries1d[15]
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL operator-(const float4x4& lhs, const float4x4& rhs)
    {
        return
        {
            lhs.entries1d[0] - rhs.entries1d[0], lhs.entries1d[1] - rhs.entries1d[1], lhs.entries1d[2] - rhs.entries1d[2], lhs.entries1d[3] - rhs.entries1d[3],
            lhs.entries1d[4] - rhs.entries1d[4], lhs.entries1d[5] - rhs.entries1d[5], lhs.entries1d[6] - rhs.entries1d[6], lhs.entries1d[7] - rhs.entries1d[7],
            lhs.entries1d[8] - rhs.entries1d[8], lhs.entries1d[9] - rhs.entries1d[9], lhs.entries1d[10] - rhs.entries1d[10], lhs.entries1d[11] - rhs.entries1d[11],
            lhs.entries1d[12] - rhs.entries1d[12], lhs.entries1d[13] - rhs.entries1d[13], lhs.entries1d[14] - rhs.entries1d[14], lhs.entries1d[15] - rhs.entries1d[15]
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL operator*(const float4x4& lhs, const float4x4& rhs)
    {
        return
        {
            (lhs.entries1d[0] * rhs.entries1d[0]) + (lhs.entries1d[1] * rhs.entries1d[4]) + (lhs.entries1d[2] * rhs.entries1d[8]) + (lhs.entries1d[3] * rhs.entries1d[12]),
            (lhs.entries1d[0] * rhs.entries1d[1]) + (lhs.entries1d[1] * rhs.entries1d[5]) + (lhs.entries1d[2] * rhs.entries1d[9]) + (lhs.entries1d[3] * rhs.entries1d[13]),
            (lhs.entries1d[0] * rhs.entries1d[2]) + (lhs.entries1d[1] * rhs.entries1d[6]) + (lhs.entries1d[2] * rhs.entries1d[10]) + (lhs.entries1d[3] * rhs.entries1d[14]),
            (lhs.entries1d[0] * rhs.entries1d[3]) + (lhs.entries1d[1] * rhs.entries1d[7]) + (lhs.entries1d[2] * rhs.entries1d[11]) + (lhs.entries1d[3] * rhs.entries1d[15]),
            (lhs.entries1d[4] * rhs.entries1d[0]) + (lhs.entries1d[5] * rhs.entries1d[4]) + (lhs.entries1d[6] * rhs.entries1d[8]) + (lhs.entries1d[7] * rhs.entries1d[12]),
            (lhs.entries1d[4] * rhs.entries1d[1]) + (lhs.entries1d[5] * rhs.entries1d[5]) + (lhs.entries1d[6] * rhs.entries1d[9]) + (lhs.entries1d[7] * rhs.entries1d[13]),
            (lhs.entries1d[4] * rhs.entries1d[2]) + (lhs.entries1d[5] * rhs.entries1d[6]) + (lhs.entries1d[6] * rhs.entries1d[10]) + (lhs.entries1d[7] * rhs.entries1d[14]),
            (lhs.entries1d[4] * rhs.entries1d[3]) + (lhs.entries1d[5] * rhs.entries1d[7]) + (lhs.entries1d[6] * rhs.entries1d[11]) + (lhs.entries1d[7] * rhs.entries1d[15]),
            (lhs.entries1d[8] * rhs.entries1d[0]) + (lhs.entries1d[9] * rhs.entries1d[4]) + (lhs.entries1d[10] * rhs.entries1d[8]) + (lhs.entries1d[11] * rhs.entries1d[12]),
            (lhs.entries1d[8] * rhs.entries1d[1]) + (lhs.entries1d[9] * rhs.entries1d[5]) + (lhs.entries1d[10] * rhs.entries1d[9]) + (lhs.entries1d[11] * rhs.entries1d[13]),
            (lhs.entries1d[8] * rhs.entries1d[2]) + (lhs.entries1d[9] * rhs.entries1d[6]) + (lhs.entries1d[10] * rhs.entries1d[10]) + (lhs.entries1d[11] * rhs.entries1d[14]),
            (lhs.entries1d[8] * rhs.entries1d[3]) + (lhs.entries1d[9] * rhs.entries1d[7]) + (lhs.entries1d[10] * rhs.entries1d[11]) + (lhs.entries1d[11] * rhs.entries1d[15]),
            (lhs.entries1d[12] * rhs.entries1d[0]) + (lhs.entries1d[13] * rhs.entries1d[4]) + (lhs.entries1d[14] * rhs.entries1d[8]) + (lhs.entries1d[15] * rhs.entries1d[12]),
            (lhs.entries1d[12] * rhs.entries1d[1]) + (lhs.entries1d[13] * rhs.entries1d[5]) + (lhs.entries1d[14] * rhs.entries1d[9]) + (lhs.entries1d[15] * rhs.entries1d[13]),
            (lhs.entries1d[12] * rhs.entries1d[2]) + (lhs.entries1d[13] * rhs.entries1d[6]) + (lhs.entries1d[14] * rhs.entries1d[10]) + (lhs.entries1d[15] * rhs.entries1d[14]),
            (lhs.entries1d[12] * rhs.entries1d[3]) + (lhs.entries1d[13] * rhs.entries1d[7]) + (lhs.entries1d[14] * rhs.entries1d[11]) + (lhs.entries1d[15] * rhs.entries1d[15])
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL operator*(const float lhs, const float4x4& rhs)
    {
        return
        {
            rhs.entries1d[ 0] * lhs, rhs.entries1d[ 1] * lhs, rhs.entries1d[ 2] * lhs, rhs.entries1d[ 3] * lhs,
            rhs.entries1d[ 4] * lhs, rhs.entries1d[ 5] * lhs, rhs.entries1d[ 6] * lhs, rhs.entries1d[ 7] * lhs,
            rhs.entries1d[ 8] * lhs, rhs.entries1d[ 9] * lhs, rhs.entries1d[10] * lhs, rhs.entries1d[11] * lhs,
            rhs.entries1d[12] * lhs, rhs.entries1d[13] * lhs, rhs.entries1d[14] * lhs, rhs.entries1d[15] * lhs
        };
    }

    GEM_INLINE float4x4 GEM_VECTORCALL operator*(const float4x4& lhs, const float rhs)
    {
        return
        {
            lhs.entries1d[0] * rhs, lhs.entries1d[1] * rhs, lhs.entries1d[2] * rhs, lhs.entries1d[3] * rhs,
            lhs.entries1d[4] * rhs, lhs.entries1d[5] * rhs, lhs.entries1d[6] * rhs, lhs.entries1d[7] * rhs,
            lhs.entries1d[8] * rhs, lhs.entries1d[9] * rhs, lhs.entries1d[10] * rhs, lhs.entries1d[11] * rhs,
            lhs.entries1d[12] * rhs, lhs.entries1d[13] * rhs, lhs.entries1d[14] * rhs, lhs.entries1d[15] * rhs
        };
    }

    GEM_INLINE float4 GEM_VECTORCALL operator*(const float4& lhs, const float4x4& rhs)
    {
        return
        {
            (lhs.components[0] * rhs.entries1d[0]) + (lhs.components[1] * rhs.entries1d[4]) + (lhs.components[2] * rhs.entries1d[8]) + (lhs.components[3] * rhs.entries1d[12]),
            (lhs.components[0] * rhs.entries1d[1]) + (lhs.components[1] * rhs.entries1d[5]) + (lhs.components[2] * rhs.entries1d[9]) + (lhs.components[3] * rhs.entries1d[13]),
            (lhs.components[0] * rhs.entries1d[2]) + (lhs.components[1] * rhs.entries1d[6]) + (lhs.components[2] * rhs.entries1d[10]) + (lhs.components[3] * rhs.entries1d[14]),
            (lhs.components[0] * rhs.entries1d[3]) + (lhs.components[1] * rhs.entries1d[7]) + (lhs.components[2] * rhs.entries1d[11]) + (lhs.components[3] * rhs.entries1d[15])
        };
    }

#pragma endregion

#pragma region double2x2

    struct double2x2
    {
        union
        {
            struct
            {
                double m00, m01;
                double m10, m11;
            };

            struct
            {
                double2 u, v;
            };

            double2 rows[2];
            double entries1d[4];
            double entries2d[2][2];
        };

        double2x2(const double u1, const double u2, const double v1, const double v2);

        double2x2(const double* p_entries);

        double2x2(const double2& u, const double2& v);

        double2x2(const double2x2& other);

        double2x2() = default;

        double2x2 transpose() const;

        double2x2 inverse() const;

        double determinant()	const;

        double2x2& GEM_VECTORCALL operator+=(const double2x2& rhs);

        double2x2& GEM_VECTORCALL operator-=(const double2x2& rhs);

        double2x2& GEM_VECTORCALL operator*=(const double rhs);

        double2x2& GEM_VECTORCALL operator=(const double2x2& rhs);

        double2& GEM_VECTORCALL operator[](const unsigned int index);

        double2 GEM_VECTORCALL operator[](const unsigned int index) const;

        double GEM_VECTORCALL operator()(const unsigned int row, const unsigned int column);

        double2x2 operator-() const;
    };

    double2x2 GEM_VECTORCALL operator+(const double2x2& lhs, const double2x2& rhs);

    double2x2 GEM_VECTORCALL operator-(const double2x2& lhs, const double2x2& rhs);

    double2x2 GEM_VECTORCALL operator*(const double2x2& lhs, const double2x2& rhs);

    double2x2 GEM_VECTORCALL operator*(const double2x2& lhs, const double rhs);

    double2x2 GEM_VECTORCALL operator*(const double lhs, const double2x2& rhs);

    double2 GEM_VECTORCALL operator*(const double2& lhs, const double2x2& rhs);

    GEM_INLINE double2x2::double2x2(const double u1, const double u2, const double v1, const double v2)
        : m00(u1), m01(u2)
        , m10(v1), m11(v2)
    {

    }

    GEM_INLINE double2x2::double2x2(const double* p_entries)
        : m00(p_entries[0]), m01(p_entries[1])
        , m10(p_entries[2]), m11(p_entries[3])
    {

    }

    GEM_INLINE double2x2::double2x2(const double2& u, const double2& v)
        : m00(u.x), m01(u.y)
        , m10(v.x), m11(v.y)
    {

    }

    GEM_INLINE double2x2::double2x2(const double2x2& o)
        : m00(o.m00), m01(o.m01)
        , m10(o.m10), m11(o.m11)
    {

    }

    GEM_INLINE double2x2 double2x2::transpose() const
    {
        return
        {
            m00, m10,
            m01, m11
        };
    }

    GEM_INLINE double2x2 double2x2::inverse() const
    {
        double2x2 out;
        double det = (entries1d[0] * entries1d[3]) - (entries1d[1] * entries1d[2]);
        double inv = 1.0f / det;
        out.entries1d[0] = +entries1d[3] * inv;
        out.entries1d[1] = -entries1d[1] * inv;
        out.entries1d[2] = -entries1d[2] * inv;
        out.entries1d[3] = +entries1d[0] * inv;
        return out;
    }

    GEM_INLINE double double2x2::determinant() const
    {
        //     0 1      
        // det 2 3 =    (0 * 3) - (1 - 2)     
        return (entries1d[0] * entries1d[3]) - (entries1d[1] * entries1d[2]);
    }

    GEM_INLINE double2x2& GEM_VECTORCALL double2x2::operator+=(const double2x2& rhs)
    {
        m00 += rhs.m00; m01 += rhs.m01;
        m10 += rhs.m10; m11 += rhs.m11;
        return *this;
    }

    GEM_INLINE double2x2& GEM_VECTORCALL double2x2::operator-=(const double2x2& rhs)
    {
        m00 -= rhs.m00; m01 -= rhs.m01;
        m10 -= rhs.m10; m11 -= rhs.m11;
        return *this;
    }

    GEM_INLINE double2x2& GEM_VECTORCALL double2x2::operator*=(const double rhs)
    {
        m00 *= rhs; m01 *= rhs;
        m10 *= rhs; m11 *= rhs;
        return *this;
    }

    GEM_INLINE double2x2& GEM_VECTORCALL double2x2::operator=(const double2x2& rhs)
    {
        m00 = rhs.m00; m01 = rhs.m01;
        m10 = rhs.m10; m11 = rhs.m11;
        return *this;
    }

    GEM_INLINE double2x2 double2x2::operator-() const
    {
        return
        {
            -m00, -m01,
            -m10, -m11
        };
    }

    GEM_INLINE double2& GEM_VECTORCALL double2x2::operator[](const unsigned int index)
    {
        return rows[index];
    }

    GEM_INLINE double2 GEM_VECTORCALL double2x2::operator[](const unsigned int index) const
    {
        return rows[index];
    }

    GEM_INLINE double GEM_VECTORCALL double2x2::operator()(const unsigned int row, const unsigned int column)
    {
        return entries2d[row][column];
    }

    GEM_INLINE double2x2 GEM_VECTORCALL operator+(const double2x2& lhs, const double2x2& rhs)
    {
        return
        {
            lhs.entries1d[0] + rhs.entries1d[0], lhs.entries1d[1] + rhs.entries1d[1],
            lhs.entries1d[2] + rhs.entries1d[2], lhs.entries1d[3] + rhs.entries1d[3]
        };
    }

    GEM_INLINE double2x2 GEM_VECTORCALL operator-(const double2x2& lhs, const double2x2& rhs)
    {
        return
        {
            lhs.entries1d[0] - rhs.entries1d[0], lhs.entries1d[1] - rhs.entries1d[1],
            lhs.entries1d[2] - rhs.entries1d[2], lhs.entries1d[3] - rhs.entries1d[3]
        };
    }

    GEM_INLINE double2x2 GEM_VECTORCALL operator*(const double2x2& lhs, const double2x2& rhs)
    {
        return
        {
            (lhs.entries1d[0] * rhs.entries1d[0]) + (lhs.entries1d[1] * rhs.entries1d[2]),
            (lhs.entries1d[0] * rhs.entries1d[1]) + (lhs.entries1d[1] * rhs.entries1d[3]),
            (lhs.entries1d[2] * rhs.entries1d[0]) + (lhs.entries1d[3] * rhs.entries1d[2]),
            (lhs.entries1d[2] * rhs.entries1d[1]) + (lhs.entries1d[3] * rhs.entries1d[3])

        };
    }

    GEM_INLINE double2 GEM_VECTORCALL operator*(const double2& lhs, const double2x2& rhs)
    {
        return
        {
            (lhs.components[0] * rhs.entries1d[0]) + (lhs.components[1] * rhs.entries1d[2]),
            (lhs.components[0] * rhs.entries1d[1]) + (lhs.components[1] * rhs.entries1d[3])
        };
    }

    GEM_INLINE double2x2 GEM_VECTORCALL operator*(const double lhs, const double2x2& rhs)
    {
        return
        {
            rhs.entries1d[0] * lhs, rhs.entries1d[1] * lhs,
            rhs.entries1d[2] * lhs, rhs.entries1d[3] * lhs
        };
    }

    GEM_INLINE double2x2 GEM_VECTORCALL operator*(const double2x2& lhs, const double rhs)
    {
        return
        {
            lhs.entries1d[0] * rhs, lhs.entries1d[1] * rhs,
            lhs.entries1d[2] * rhs, lhs.entries1d[3] * rhs
        };
    }

#pragma endregion

#pragma region double3x3

    struct double3x3
    {
        union
        {
            struct
            {
                double m00, m01, m02;
                double m10, m11, m12;
                double m20, m21, m22;
            };

            struct
            {
                double3 u, v, w;
            };

            double3 rows[3];
            double entries1d[9];
            double entries2d[3][3];
        };

        static double3x3 identity();

        static double3x3 GEM_VECTORCALL rotate_euler(const double pitch, const double yaw, const double roll);

        static double3x3 GEM_VECTORCALL rotate_axis_angle(const double3& a, const double radians);

        static double3x3 GEM_VECTORCALL rotate_axis_x(const double radians);

        static double3x3 GEM_VECTORCALL rotate_axis_y(const double radians);

        static double3x3 GEM_VECTORCALL rotate_axis_z(const double radians);

        static double3x3 GEM_VECTORCALL rotate_forward_lh(const double3& direction, const double3& up);

        static double3x3 GEM_VECTORCALL rotate_forward_rh(const double3& direction, const double3& up);

        static double3x3 GEM_VECTORCALL scale(const double3& s);

        double3x3(const double u1, const double u2, const double u3, const double v1, const double v2, const double v3, const double w1, const double w2, const double w3);

        double3x3(const double* p_entries);

        double3x3(const double3& u, const double3& v, const double3& w);

        double3x3(const double3x3& o);

        double3x3() = default;

        double3x3 transpose() const;

        double3x3 inverse() const;

        double determinant()	const;

        double3 euler() const;

        double3x3& GEM_VECTORCALL operator=(const double3x3& rhs);

        double3x3& GEM_VECTORCALL operator+=(const double3x3& rhs);

        double3x3& GEM_VECTORCALL operator-=(const double3x3& rhs);

        double3x3& GEM_VECTORCALL operator*=(const double rhs);

        double3& GEM_VECTORCALL operator[](const unsigned int index);

        double3 GEM_VECTORCALL operator[](const unsigned int index) const;

        double GEM_VECTORCALL operator()(const unsigned int row, const unsigned int column);

        double3x3 operator-() const;
    };

    double3x3 GEM_VECTORCALL operator+(const double3x3& lhs, const double3x3& rhs);

    double3x3 GEM_VECTORCALL operator-(const double3x3& lhs, const double3x3& rhs);

    double3x3 GEM_VECTORCALL operator*(const double3x3& lhs, const double3x3& rhs);

    double3x3 GEM_VECTORCALL operator*(const double3x3& lhs, const double rhs);

    double3x3 GEM_VECTORCALL operator*(const double lhs, const double3x3& rhs);

    double3 GEM_VECTORCALL operator*(const double3& lhs, const double3x3& rhs);

    GEM_INLINE double3x3 double3x3::identity()
    {
        return
        {
            1.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 1.0f
        };
    }

    GEM_INLINE double3x3 GEM_VECTORCALL double3x3::rotate_euler(const double pitch, const double yaw, const double roll)
    {
        double ch = cosf(yaw); double cp = cosf(pitch); double cb = cosf(roll);
        double sh = sinf(yaw); double sp = sinf(pitch); double sb = sinf(roll);

        return
        {
             ch * cb + sh * sp * sb, sb * cp, -sh * cb + ch * sp * sb,
            -ch * sb + sh * sp * cb, cb * cp,  sb * sh + ch * sp * cb,
             sh * cp, -sp,  ch * cp
        };
    }

    GEM_INLINE double3x3 GEM_VECTORCALL double3x3::rotate_axis_angle(const double3& a, const double radians)
    {
        double s = sinf(radians);
        double c = cosf(radians);

        return
        {
            c + ((1 - c) * a.x * a.x), ((1 - c) * a.x * a.y) + (s * a.z), ((1 - c) * a.x * a.z) - (s * a.y),
            ((1 - c) * a.x * a.y) - (s * a.z), c + ((1 - c) * a.y * a.y), ((1 - c) * a.y * a.z) + (s * a.x),
            ((1 - c) * a.x * a.z) + (s * a.y), ((1 - c) * a.y * a.z) - (s * a.x), c + ((1 - c) * a.z * a.z)
        };
    }

    GEM_INLINE double3x3 GEM_VECTORCALL double3x3::rotate_axis_x(const double radians)
    {
        double s = sinf(radians);
        double c = cosf(radians);
        return
        {
            1.0f, 0.0f, 0.0f,
            0.0f, c, s,
            0.0f, -s, c
        };
    }

    GEM_INLINE double3x3 GEM_VECTORCALL double3x3::rotate_axis_y(const double radians)
    {
        double s = sin(radians);
        double c = cos(radians);
        return
        {
            c, 0.0f, -s,
            0.0f, 1.0f, 0.0f,
            s, 0.0f, c
        };
    }

    GEM_INLINE double3x3 GEM_VECTORCALL double3x3::rotate_axis_z(const double radians)
    {
        double s = sin(radians);
        double c = cos(radians);
        return
        {
            c, s, 0.0f,
            -s, c, 0.0f,
            0.0f, 0.0f, 1.0f
        };
    }

    GEM_INLINE double3x3 GEM_VECTORCALL double3x3::rotate_forward_lh(const double3& direction, const double3& up)
    {
        double3 f = normalize(direction);
        double3 s = normalize(cross(up, f));
        double3 u = cross(f, s);

        return
        {
            s.x, s.y, s.z,
            u.x, u.y, u.z,
            f.x, f.y, f.z
        };
    }

    GEM_INLINE double3x3 GEM_VECTORCALL double3x3::rotate_forward_rh(const double3& direction, const double3& up)
    {
        return rotate_forward_lh(-direction, up);
    }

    GEM_INLINE double3x3 GEM_VECTORCALL double3x3::scale(const double3& s)
    {
        return
        {
            s.x, 0.0f, 0.0f,
            0.0f, s.y, 0.0f,
            0.0f, 0.0f, s.z
        };
    }

    GEM_INLINE double3x3::double3x3(
        const double u1, const double u2, const double u3,
        const double v1, const double v2, const double v3,
        const double w1, const double w2, const double w3)
        : m00(u1), m01(u2), m02(u3)
        , m10(v1), m11(v2), m12(v3)
        , m20(w1), m21(w2), m22(w3)
    {

    }

    GEM_INLINE double3x3::double3x3(const double* o)
        : m00(o[0]), m01(o[1]), m02(o[2])
        , m10(o[3]), m11(o[4]), m12(o[5])
        , m20(o[6]), m21(o[7]), m22(o[8])
    {

    }

    GEM_INLINE double3x3::double3x3(const double3& u, const double3& v, const double3& w)
        : m00(u.x), m01(u.y), m02(u.z)
        , m10(v.x), m11(v.y), m12(v.z)
        , m20(w.x), m21(w.y), m22(w.z)
    {

    };

    GEM_INLINE double3x3::double3x3(const double3x3& o)
        : m00(o.m00), m01(o.m01), m02(o.m02)
        , m10(o.m10), m11(o.m11), m12(o.m12)
        , m20(o.m20), m21(o.m21), m22(o.m22)
    {

    }

    GEM_INLINE double3x3 double3x3::transpose() const
    {
        return
        {
            m00, m10, m20,
            m01, m11, m21,
            m02, m12, m22
        };
    }

    GEM_INLINE double3x3 double3x3::inverse() const
    {
        double3x3 out;

        double det =
            (((entries1d[1] * entries1d[5]) - (entries1d[2] * entries1d[4])) * entries1d[6]) +
            (((entries1d[2] * entries1d[3]) - (entries1d[0] * entries1d[5])) * entries1d[7]) +
            (((entries1d[0] * entries1d[4]) - (entries1d[1] * entries1d[3])) * entries1d[8]);

        double c11 = +((entries1d[4] * entries1d[8]) - (entries1d[5] * entries1d[7]));
        double c12 = -((entries1d[3] * entries1d[8]) - (entries1d[5] * entries1d[6]));
        double c13 = +((entries1d[3] * entries1d[7]) - (entries1d[4] * entries1d[6]));
        double c21 = -((entries1d[1] * entries1d[8]) - (entries1d[2] * entries1d[7]));
        double c22 = +((entries1d[0] * entries1d[8]) - (entries1d[2] * entries1d[6]));
        double c23 = -((entries1d[0] * entries1d[7]) - (entries1d[1] * entries1d[6]));
        double c31 = +((entries1d[1] * entries1d[5]) - (entries1d[2] * entries1d[4]));
        double c32 = -((entries1d[0] * entries1d[5]) - (entries1d[2] * entries1d[3]));
        double c33 = +((entries1d[0] * entries1d[4]) - (entries1d[1] * entries1d[3]));

        double inv = 1.0f / det;
        out.entries1d[0] = c11 * inv;
        out.entries1d[1] = c21 * inv;
        out.entries1d[2] = c31 * inv;
        out.entries1d[3] = c12 * inv;
        out.entries1d[4] = c22 * inv;
        out.entries1d[5] = c32 * inv;
        out.entries1d[6] = c13 * inv;
        out.entries1d[7] = c23 * inv;
        out.entries1d[8] = c33 * inv;
        return out;
    }

    GEM_INLINE double double3x3::determinant() const
    {
        // ((row0 x row1) * row2)
        // (ay*bz - az*by) * cx
        // (az*bx - ax*bz) * cy
        // (ax*by - ay*bx) * cz

        return
            (((entries1d[1] * entries1d[5]) - (entries1d[2] * entries1d[4])) * entries1d[6]) +
            (((entries1d[2] * entries1d[3]) - (entries1d[0] * entries1d[5])) * entries1d[7]) +
            (((entries1d[0] * entries1d[4]) - (entries1d[1] * entries1d[3])) * entries1d[8]);
    }

    GEM_INLINE double3 double3x3::euler() const
    {
        double3 out;
        double sp = -entries1d[7];
        if (sp <= 1.0f)
        {
            out.x = -1.570796f;
        }
        else if (sp >= 1.0f)
        {
            out.x = -1.570796f;
        }
        else
        {
            out.x = asin(sp);
        }

        if (abs(sp) > 0.9999f)
        {
            out.z = 0.0f;
            out.y = atan2(-entries1d[2], entries1d[0]);
        }
        else
        {
            out.y = atan2(entries1d[6], entries1d[8]);
            out.z = atan2(entries1d[1], entries1d[4]);
        }
        return out;
    }

    GEM_INLINE double3x3& GEM_VECTORCALL double3x3::operator=(const double3x3& rhs)
    {
        m00 = rhs.m00; m01 = rhs.m01; m02 = rhs.m02;
        m10 = rhs.m10; m11 = rhs.m11; m12 = rhs.m12;
        m20 = rhs.m20; m21 = rhs.m21; m22 = rhs.m22;
        return *this;
    }

    GEM_INLINE double3x3& GEM_VECTORCALL double3x3::operator+=(const double3x3& rhs)
    {
        m00 += rhs.m00; m01 += rhs.m01; m02 += rhs.m02;
        m10 += rhs.m10; m11 += rhs.m11; m12 += rhs.m12;
        m20 += rhs.m20; m21 += rhs.m21; m22 += rhs.m22;
        return *this;
    }

    GEM_INLINE double3x3& GEM_VECTORCALL double3x3::operator-=(const double3x3& rhs)
    {
        m00 -= rhs.m00; m01 -= rhs.m01; m02 -= rhs.m02;
        m10 -= rhs.m10; m11 -= rhs.m11; m12 -= rhs.m12;
        m20 -= rhs.m20; m21 -= rhs.m21; m22 -= rhs.m22;
        return *this;
    }

    GEM_INLINE double3x3& GEM_VECTORCALL double3x3::operator*=(const double rhs)
    {
        m00 *= rhs; m01 *= rhs; m02 *= rhs;
        m10 *= rhs; m11 *= rhs; m12 *= rhs;
        m20 *= rhs; m21 *= rhs; m22 *= rhs;
        return *this;
    }

    GEM_INLINE double3& GEM_VECTORCALL double3x3::operator[](const unsigned int index)
    {
        return rows[index];
    }

    GEM_INLINE double3 GEM_VECTORCALL double3x3::operator[](const unsigned int index) const
    {
        return rows[index];
    }

    GEM_INLINE double GEM_VECTORCALL double3x3::operator()(const unsigned int row, const unsigned int column)
    {
        return entries2d[row][column];
    }

    GEM_INLINE double3x3 double3x3::operator-() const
    {
        return
        {
            -m00, -m01, -m02,
            -m10, -m11, -m12,
            -m20, -m21, -m22
        };
    }

    GEM_INLINE double3x3 GEM_VECTORCALL operator+(const double3x3& lhs, const double3x3& rhs)
    {
        return
        {
            lhs.entries1d[0] + rhs.entries1d[0], lhs.entries1d[1] + rhs.entries1d[1], lhs.entries1d[2] + rhs.entries1d[2],
            lhs.entries1d[3] + rhs.entries1d[3], lhs.entries1d[4] + rhs.entries1d[4], lhs.entries1d[5] + rhs.entries1d[5],
            lhs.entries1d[6] + rhs.entries1d[6], lhs.entries1d[7] + rhs.entries1d[7], lhs.entries1d[8] + rhs.entries1d[8],
        };
    }

    GEM_INLINE double3x3 GEM_VECTORCALL operator-(const double3x3& lhs, const double3x3& rhs)
    {
        return
        {
            lhs.entries1d[0] - rhs.entries1d[0], lhs.entries1d[1] - rhs.entries1d[1], lhs.entries1d[2] - rhs.entries1d[2],
            lhs.entries1d[3] - rhs.entries1d[3], lhs.entries1d[4] - rhs.entries1d[4], lhs.entries1d[5] - rhs.entries1d[5],
            lhs.entries1d[6] - rhs.entries1d[6], lhs.entries1d[7] - rhs.entries1d[7], lhs.entries1d[8] - rhs.entries1d[8],
        };
    }

    GEM_INLINE double3x3 GEM_VECTORCALL operator*(const double3x3& lhs, const double3x3& rhs)
    {
        return
        {
            (lhs.entries1d[0] * rhs.entries1d[0]) + (lhs.entries1d[1] * rhs.entries1d[3]) + (lhs.entries1d[2] * rhs.entries1d[6]),
            (lhs.entries1d[0] * rhs.entries1d[1]) + (lhs.entries1d[1] * rhs.entries1d[4]) + (lhs.entries1d[2] * rhs.entries1d[7]),
            (lhs.entries1d[0] * rhs.entries1d[2]) + (lhs.entries1d[1] * rhs.entries1d[5]) + (lhs.entries1d[2] * rhs.entries1d[8]),
            (lhs.entries1d[3] * rhs.entries1d[0]) + (lhs.entries1d[4] * rhs.entries1d[3]) + (lhs.entries1d[5] * rhs.entries1d[6]),
            (lhs.entries1d[3] * rhs.entries1d[1]) + (lhs.entries1d[4] * rhs.entries1d[4]) + (lhs.entries1d[5] * rhs.entries1d[7]),
            (lhs.entries1d[3] * rhs.entries1d[2]) + (lhs.entries1d[4] * rhs.entries1d[5]) + (lhs.entries1d[5] * rhs.entries1d[8]),
            (lhs.entries1d[6] * rhs.entries1d[0]) + (lhs.entries1d[7] * rhs.entries1d[3]) + (lhs.entries1d[8] * rhs.entries1d[6]),
            (lhs.entries1d[6] * rhs.entries1d[1]) + (lhs.entries1d[7] * rhs.entries1d[4]) + (lhs.entries1d[8] * rhs.entries1d[7]),
            (lhs.entries1d[6] * rhs.entries1d[2]) + (lhs.entries1d[7] * rhs.entries1d[5]) + (lhs.entries1d[8] * rhs.entries1d[8])
        };
    }

    GEM_INLINE double3 GEM_VECTORCALL operator*(const double3& lhs, const double3x3& rhs)
    {
        return
        {

            (lhs.components[0] * rhs.entries1d[0]) + (lhs.components[1] * rhs.entries1d[3]) + (lhs.components[2] * rhs.entries1d[6]),
            (lhs.components[0] * rhs.entries1d[1]) + (lhs.components[1] * rhs.entries1d[4]) + (lhs.components[2] * rhs.entries1d[7]),
            (lhs.components[0] * rhs.entries1d[2]) + (lhs.components[1] * rhs.entries1d[5]) + (lhs.components[2] * rhs.entries1d[8])
        };
    }

    GEM_INLINE double3x3 GEM_VECTORCALL operator*(const double lhs, const double3x3& rhs)
    {
        return
        {
            rhs.entries1d[0] * lhs, rhs.entries1d[1] * lhs, rhs.entries1d[2] * lhs,
            rhs.entries1d[3] * lhs, rhs.entries1d[4] * lhs, rhs.entries1d[5] * lhs,
            rhs.entries1d[6] * lhs, rhs.entries1d[7] * lhs, rhs.entries1d[8] * lhs
        };
    }

    GEM_INLINE double3x3 GEM_VECTORCALL operator*(const double3x3& lhs, const double rhs)
    {
        return
        {
            lhs.entries1d[0] * rhs, lhs.entries1d[1] * rhs, lhs.entries1d[2] * rhs,
            lhs.entries1d[3] * rhs, lhs.entries1d[4] * rhs, lhs.entries1d[5] * rhs,
            lhs.entries1d[6] * rhs, lhs.entries1d[7] * rhs, lhs.entries1d[8] * rhs
        };
    }

#pragma endregion

#pragma region double4x3

    struct double4x3
    {
        union
        {
            struct
            {
                double m00, m01, m02;
                double m10, m11, m12;
                double m20, m21, m22;
                double m30, m31, m32;
            };

            struct
            {
                double3 u, v, w, t;
            };

            double3 rows[4];
            double entries1d[12];
            double entries2d[4][3];
        };

        static double4x3 identity();

        static double4x3 GEM_VECTORCALL rotate_euler(const double pitch, const double yaw, const double roll);

        static double4x3 GEM_VECTORCALL rotate_axis_angle(const double3& a, const double angle);

        static double4x3 GEM_VECTORCALL rotate_axis_x(const double angle);

        static double4x3 GEM_VECTORCALL rotate_axis_y(const double angle);

        static double4x3 GEM_VECTORCALL rotate_axis_z(const double angle);

        static double4x3 GEM_VECTORCALL rotate_forward_lh(const double3& direction, const double3& position, const double3& up);

        static double4x3 GEM_VECTORCALL rotate_forward_rh(const double3& direction, const double3& position, const double3& up);

        static double4x3 GEM_VECTORCALL scale(const double3& s);

        static double4x3 GEM_VECTORCALL translate(const double3& t);

        double4x3(const double u1, const double u2, const double u3,
            const double v1, const double v2, const double v3,
            const double w1, const double w2, const double w3,
            const double t1, const double t2, const double t3);

        double4x3(const double* o);

        double4x3(const double3& u,
            const double3& v,
            const double3& w,
            const double3& t);

        double4x3(const double4x3& o);

        double4x3() = default;

        double4x3 transpose() const;

        double4x3 inverse() const;

        double determinant()	const;

        double3 euler() const;

        double4x3& GEM_VECTORCALL operator=(const double4x3& rhs);

        double4x3 GEM_VECTORCALL operator+=(const double4x3& rhs);

        double4x3 GEM_VECTORCALL operator-=(const double4x3& rhs);

        double4x3 GEM_VECTORCALL operator*=(const double rhs);

        double3& GEM_VECTORCALL operator[](const unsigned int index);

        double3 GEM_VECTORCALL operator[](const unsigned int index) const;

        double GEM_VECTORCALL operator()(const unsigned int row, const unsigned int column);

        double4x3 operator-() const;
    };

    double4x3 GEM_VECTORCALL operator+(const double4x3& lhs, const double4x3& rhs);

    double4x3 GEM_VECTORCALL operator-(const double4x3& lhs, const double4x3& rhs);

    double4x3 GEM_VECTORCALL operator*(const double4x3& lhs, const double4x3& rhs);

    double4x3 GEM_VECTORCALL operator*(const double4x3& lhs, const double rhs);

    double4x3 GEM_VECTORCALL operator*(const double lhs, const double4x3& rhs);

    double3 GEM_VECTORCALL operator*(const double3& lhs, const double4x3& rhs);

    GEM_INLINE double4x3 double4x3::identity()
    {
        return
        {
            1.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 1.0f,
            0.0f, 0.0f, 0.0f
        };
    }

    GEM_INLINE double4x3 GEM_VECTORCALL double4x3::rotate_euler(const double pitch, const double yaw, const double roll)
    {
        double ch = cosf(yaw); double cp = cosf(pitch); double cb = cosf(roll);
        double sh = sinf(yaw); double sp = sinf(pitch); double sb = sinf(roll);

        return
        {
             ch * cb + sh * sp * sb, sb * cp, -sh * cb + ch * sp * sb,
            -ch * sb + sh * sp * cb, cb * cp,  sb * sh + ch * sp * cb,
             sh * cp, -sp,  ch * cp,
             0.0f, 0.0f, 0.0f
        };
    }

    GEM_INLINE double4x3 GEM_VECTORCALL double4x3::rotate_axis_angle(const double3& a, const double angle)
    {
        double s = sin(angle);
        double c = cos(angle);
        double ax2 = a.x * a.x;

        return
        {
            c + ((1 - c) * a.x * a.x), ((1 - c) * a.x * a.y) + (s * a.z), ((1 - c) * a.x * a.z) - (s * a.y),
            ((1 - c) * a.x * a.y) - (s * a.z), c + ((1 - c) * a.y * a.y), ((1 - c) * a.y * a.z) + (s * a.x),
            ((1 - c) * a.x * a.z) + (s * a.y), ((1 - c) * a.y * a.z) - (s * a.x), c + ((1 - c) * a.z * a.z),
            0.0f, 0.0f, 0.0
        };
    }

    GEM_INLINE double4x3 GEM_VECTORCALL double4x3::rotate_axis_x(const double angle)
    {
        double s = sin(angle);
        double c = cos(angle);
        return
        {
            1.0f, 0.0f, 0.0f,
            0.0f, c, s,
            0.0f, -s, c,
            0.0f, 0.0f, 0.0
        };
    }

    GEM_INLINE  double4x3 GEM_VECTORCALL double4x3::rotate_axis_y(const double angle)
    {
        double s = sin(angle);
        double c = cos(angle);
        return
        {
            c, 0.0f, -s,
            0.0f, 1.0f, 0.0f,
            s, 0.0f, c,
            0.0f, 0.0f, 0.0
        };
    }

    GEM_INLINE double4x3 GEM_VECTORCALL double4x3::rotate_axis_z(const double angle)
    {
        double s = sin(angle);
        double c = cos(angle);
        return
        {
            c, s, 0.0f,
            -s, c, 0.0f,
            0.0f, 0.0f, 1.0f,
            0.0f, 0.0f, 0.0
        };
    }

    GEM_INLINE double4x3 GEM_VECTORCALL double4x3::rotate_forward_lh(const double3& direction, const double3& position, const double3& up)
    {
        double3 f = normalize(direction);
        double3 s = normalize(cross(up, f));
        double3 u = cross(f, s);

        return
        {
            s.x, s.y, s.z,
            u.x, u.y, u.z,
            f.x, f.y, f.z,
            position.x, position.y, position.z
        };
    }

    GEM_INLINE double4x3 GEM_VECTORCALL double4x3::rotate_forward_rh(const double3& direction, const double3& position, const double3& up)
    {
        return rotate_forward_lh(-direction, position, up);
    }

    GEM_INLINE double4x3 GEM_VECTORCALL double4x3::scale(const double3& s)
    {
        return
        {
            s.x, 0.0f, 0.0f,
            0.0f, s.y, 0.0f,
            0.0f, 0.0f, s.z,
            0.0f, 0.0f, 0.0
        };
    }

    GEM_INLINE double4x3 GEM_VECTORCALL double4x3::translate(const double3& t)
    {
        return
        {
            1.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 1.0f,
            t.x, t.y, t.z
        };
    }

    GEM_INLINE double4x3::double4x3(
        const double u1, const double u2, const double u3,
        const double v1, const double v2, const double v3,
        const double w1, const double w2, const double w3,
        const double t1, const double t2, const double t3)
        : m00(u1), m01(u2), m02(u3)
        , m10(v1), m11(v2), m12(v3)
        , m20(w1), m21(w2), m22(w3)
        , m30(t1), m31(t2), m32(t3)
    {

    }

    GEM_INLINE double4x3::double4x3(const double* o)
        : m00(o[0]), m01(o[1]), m02(o[2])
        , m10(o[3]), m11(o[4]), m12(o[5])
        , m20(o[6]), m21(o[7]), m22(o[8])
        , m30(o[9]), m31(o[10]), m32(o[11])
    {

    }

    GEM_INLINE double4x3::double4x3(const double3& u, const double3& v, const double3& w, const double3& t)
        : m00(u.x), m01(u.y), m02(u.z)
        , m10(v.x), m11(v.y), m12(v.z)
        , m20(w.x), m21(w.y), m22(w.z)
        , m30(t.x), m31(t.y), m32(t.z)
    {

    };

    GEM_INLINE double4x3::double4x3(const double4x3& o)
        : m00(o.m00), m01(o.m01), m02(o.m02)
        , m10(o.m10), m11(o.m11), m12(o.m12)
        , m20(o.m20), m21(o.m21), m22(o.m22)
        , m30(o.m30), m31(o.m31), m32(o.m32)
    {

    }

    GEM_INLINE double4x3 double4x3::transpose() const
    {
        return
        {
            m00, m10, m20, m30,
            m01, m11, m12, m31,
            m02, m12, m22, m32
        };
    }

    GEM_INLINE double4x3 double4x3::inverse() const
    {
        double4x3 out;
        double c11 = +((entries1d[4] * entries1d[8]) - (entries1d[5] * entries1d[7]));
        double c12 = -((entries1d[3] * entries1d[8]) - (entries1d[5] * entries1d[6]));
        double c13 = +((entries1d[3] * entries1d[7]) - (entries1d[4] * entries1d[6]));
        double c14 = -((((entries1d[4] * entries1d[8]) - (entries1d[5] * entries1d[7])) * entries1d[9]) + (((entries1d[5] * entries1d[6]) - (entries1d[3] * entries1d[8])) * entries1d[10]) + (((entries1d[3] * entries1d[7]) - (entries1d[4] * entries1d[6])) * entries1d[11]));

        double det = (entries1d[0] * c11) + (entries1d[1] * c12) + (entries1d[2] * c13);

        double c21 = -((entries1d[1] * entries1d[8]) - (entries1d[2] * entries1d[7]));
        double c22 = +((entries1d[0] * entries1d[8]) - (entries1d[2] * entries1d[6]));
        double c23 = -((entries1d[0] * entries1d[7]) - (entries1d[1] * entries1d[6]));
        double c24 = +((((entries1d[1] * entries1d[8]) - (entries1d[2] * entries1d[7])) * entries1d[9]) + (((entries1d[2] * entries1d[6]) - (entries1d[0] * entries1d[8])) * entries1d[10]) + (((entries1d[0] * entries1d[7]) - (entries1d[1] * entries1d[6])) * entries1d[11]));

        double c31 = +((entries1d[1] * entries1d[5]) - (entries1d[2] * entries1d[4]));
        double c32 = -((entries1d[0] * entries1d[5]) - (entries1d[2] * entries1d[3]));
        double c33 = +((entries1d[0] * entries1d[4]) - (entries1d[1] * entries1d[3]));
        double c34 = -((((entries1d[1] * entries1d[5]) - (entries1d[2] * entries1d[4])) * entries1d[9]) + (((entries1d[2] * entries1d[3]) - (entries1d[0] * entries1d[5])) * entries1d[10]) + (((entries1d[0] * entries1d[4]) - (entries1d[1] * entries1d[3])) * entries1d[11]));

        // double c41 = double(0);
        // double c42 = double(0);
        // double c43 = double(0);
        double c44 = +((((entries1d[1] * entries1d[5]) - (entries1d[2] * entries1d[4])) * entries1d[6]) + (((entries1d[2] * entries1d[3]) - (entries1d[0] * entries1d[5])) * entries1d[7]) + (((entries1d[0] * entries1d[4]) - (entries1d[1] * entries1d[3])) * entries1d[8]));

        double inv = 1.0f / det;

        out.entries1d[0] = c11 * inv;
        out.entries1d[1] = c21 * inv;
        out.entries1d[2] = c31 * inv;
        out.entries1d[3] = c12 * inv;
        out.entries1d[4] = c22 * inv;
        out.entries1d[5] = c32 * inv;
        out.entries1d[6] = c13 * inv;
        out.entries1d[7] = c23 * inv;
        out.entries1d[8] = c33 * inv;
        out.entries1d[9] = c14 * inv;
        out.entries1d[10] = c24 * inv;
        out.entries1d[11] = c34 * inv;
        return out;
    }

    GEM_INLINE double double4x3::determinant() const
    {
        double c11 = +((entries1d[4] * entries1d[8]) - (entries1d[5] * entries1d[7]));
        double c12 = -((entries1d[3] * entries1d[8]) - (entries1d[5] * entries1d[6]));
        double c13 = +((entries1d[3] * entries1d[7]) - (entries1d[4] * entries1d[6]));
        double c14 = -((((entries1d[4] * entries1d[8]) - (entries1d[5] * entries1d[7])) * entries1d[9]) + (((entries1d[5] * entries1d[6]) - (entries1d[3] * entries1d[8])) * entries1d[10]) + (((entries1d[3] * entries1d[7]) - (entries1d[4] * entries1d[6])) * entries1d[11]));

        return (entries1d[0] * c11) + (entries1d[1] * c12) + (entries1d[2] * c13);
    }

    GEM_INLINE double3 double4x3::euler() const
    {
        double3 out;
        double sp = -entries1d[7];
        if (sp <= 1.0f)
        {
            out.x = -1.570796f;
        }
        else if (sp >= 1.0f)
        {
            out.x = -1.570796f;
        }
        else
        {
            out.x = asin(sp);
        }

        if (abs(sp) > 0.9999f)
        {
            out.z = 0.0f;
            out.y = atan2(-entries1d[2], entries1d[0]);
        }
        else
        {
            out.y = atan2(entries1d[6], entries1d[8]);
            out.z = atan2(entries1d[1], entries1d[4]);
        }
        return out;
    }

    GEM_INLINE double4x3& GEM_VECTORCALL double4x3::operator=(const double4x3& rhs)
    {
        m00 = rhs.m00; m01 = rhs.m01; m02 = rhs.m02;
        m10 = rhs.m10; m11 = rhs.m11; m12 = rhs.m12;
        m20 = rhs.m20; m21 = rhs.m21; m22 = rhs.m22;
        m30 = rhs.m30; m31 = rhs.m31; m32 = rhs.m32;
        return *this;
    }

    GEM_INLINE double4x3 GEM_VECTORCALL double4x3::operator+=(const double4x3& rhs)
    {
        m00 += rhs.m00; m01 += rhs.m01; m02 += rhs.m02;
        m10 += rhs.m10; m11 += rhs.m11; m12 += rhs.m12;
        m20 += rhs.m20; m21 += rhs.m21; m22 += rhs.m22;
        m30 += rhs.m30; m31 += rhs.m31; m32 += rhs.m32;
        return *this;
    }

    GEM_INLINE double4x3 GEM_VECTORCALL double4x3::operator-=(const double4x3& rhs)
    {
        m00 -= rhs.m00; m01 -= rhs.m01; m02 -= rhs.m02;
        m10 -= rhs.m10; m11 -= rhs.m11; m12 -= rhs.m12;
        m20 -= rhs.m20; m21 -= rhs.m21; m22 -= rhs.m22;
        m30 -= rhs.m30; m31 -= rhs.m31; m32 -= rhs.m32;
        return *this;
    }

    GEM_INLINE double4x3 GEM_VECTORCALL double4x3::operator*=(const double rhs)
    {
        m00 *= rhs; m01 *= rhs; m02 *= rhs;
        m10 *= rhs; m11 *= rhs; m12 *= rhs;
        m20 *= rhs; m21 *= rhs; m22 *= rhs;
        m30 *= rhs; m31 *= rhs; m32 *= rhs;
        return *this;
    }

    GEM_INLINE double3& GEM_VECTORCALL double4x3::operator[](const unsigned int index)
    {
        return rows[index];
    }

    GEM_INLINE double3 GEM_VECTORCALL double4x3::operator[](const unsigned int index) const
    {
        return rows[index];
    }

    GEM_INLINE double GEM_VECTORCALL double4x3::operator()(const unsigned int row, const unsigned int column)
    {
        return entries2d[row][column];
    }

    GEM_INLINE double4x3 double4x3::operator-() const
    {
        return
        {
            -m00, -m01, -m02,
            -m10, -m11, -m12,
            -m20, -m21, -m22,
            -m30, -m31, -m32
        };
    }

    GEM_INLINE double4x3 GEM_VECTORCALL operator+(const double4x3& lhs, const double4x3& rhs)
    {
        return
        {
            lhs.entries1d[0] + rhs.entries1d[0], lhs.entries1d[1] + rhs.entries1d[1], lhs.entries1d[2] + rhs.entries1d[2],
            lhs.entries1d[3] + rhs.entries1d[3], lhs.entries1d[4] + rhs.entries1d[4], lhs.entries1d[5] + rhs.entries1d[5],
            lhs.entries1d[6] + rhs.entries1d[6], lhs.entries1d[7] + rhs.entries1d[7], lhs.entries1d[8] + rhs.entries1d[8],
            lhs.entries1d[9] + rhs.entries1d[9], lhs.entries1d[10] + rhs.entries1d[10], lhs.entries1d[11] + rhs.entries1d[11]
        };
    }

    GEM_INLINE double4x3 GEM_VECTORCALL operator-(const double4x3& lhs, const double4x3& rhs)
    {
        return
        {
            lhs.entries1d[0] - rhs.entries1d[0], lhs.entries1d[1] - rhs.entries1d[1], lhs.entries1d[2] - rhs.entries1d[2],
            lhs.entries1d[3] - rhs.entries1d[3], lhs.entries1d[4] - rhs.entries1d[4], lhs.entries1d[5] - rhs.entries1d[5],
            lhs.entries1d[6] - rhs.entries1d[6], lhs.entries1d[7] - rhs.entries1d[7], lhs.entries1d[8] - rhs.entries1d[8],
            lhs.entries1d[9] - rhs.entries1d[9], lhs.entries1d[10] - rhs.entries1d[10], lhs.entries1d[11] - rhs.entries1d[11]
        };
    }

    GEM_INLINE double4x3 GEM_VECTORCALL operator*(const double4x3& lhs, const double4x3& rhs)
    {
        return
        {
            (lhs.entries1d[0] * rhs.entries1d[0]) + (lhs.entries1d[1] * rhs.entries1d[3]) + (lhs.entries1d[2] * rhs.entries1d[6]),
            (lhs.entries1d[0] * rhs.entries1d[1]) + (lhs.entries1d[1] * rhs.entries1d[4]) + (lhs.entries1d[2] * rhs.entries1d[7]),
            (lhs.entries1d[0] * rhs.entries1d[2]) + (lhs.entries1d[1] * rhs.entries1d[5]) + (lhs.entries1d[2] * rhs.entries1d[8]),
            (lhs.entries1d[3] * rhs.entries1d[0]) + (lhs.entries1d[4] * rhs.entries1d[3]) + (lhs.entries1d[5] * rhs.entries1d[6]),
            (lhs.entries1d[3] * rhs.entries1d[1]) + (lhs.entries1d[4] * rhs.entries1d[4]) + (lhs.entries1d[5] * rhs.entries1d[7]),
            (lhs.entries1d[3] * rhs.entries1d[2]) + (lhs.entries1d[4] * rhs.entries1d[5]) + (lhs.entries1d[5] * rhs.entries1d[8]),
            (lhs.entries1d[6] * rhs.entries1d[0]) + (lhs.entries1d[7] * rhs.entries1d[3]) + (lhs.entries1d[8] * rhs.entries1d[6]),
            (lhs.entries1d[6] * rhs.entries1d[1]) + (lhs.entries1d[7] * rhs.entries1d[4]) + (lhs.entries1d[8] * rhs.entries1d[7]),
            (lhs.entries1d[6] * rhs.entries1d[2]) + (lhs.entries1d[7] * rhs.entries1d[5]) + (lhs.entries1d[8] * rhs.entries1d[8]),
            (lhs.entries1d[9] * rhs.entries1d[0]) + (lhs.entries1d[10] * rhs.entries1d[3]) + (lhs.entries1d[11] * rhs.entries1d[6]) + rhs.entries1d[9],
            (lhs.entries1d[9] * rhs.entries1d[1]) + (lhs.entries1d[10] * rhs.entries1d[4]) + (lhs.entries1d[11] * rhs.entries1d[7]) + rhs.entries1d[10],
            (lhs.entries1d[9] * rhs.entries1d[2]) + (lhs.entries1d[10] * rhs.entries1d[5]) + (lhs.entries1d[11] * rhs.entries1d[8]) + rhs.entries1d[11]
        };
    }

    GEM_INLINE double4x3 GEM_VECTORCALL operator*(const double lhs, const double4x3& rhs)
    {
        return
        {
            rhs.entries1d[0] * lhs, rhs.entries1d[1] * lhs, rhs.entries1d[2] * lhs,
            rhs.entries1d[3] * lhs, rhs.entries1d[4] * lhs, rhs.entries1d[5] * lhs,
            rhs.entries1d[6] * lhs, rhs.entries1d[7] * lhs, rhs.entries1d[8] * lhs,
            rhs.entries1d[9] * lhs, rhs.entries1d[10] * lhs, rhs.entries1d[11] * lhs
        };
    }

    GEM_INLINE double4x3 GEM_VECTORCALL operator*(const double4x3& lhs, const double rhs)
    {
        return
        {
            lhs.entries1d[0] * rhs, lhs.entries1d[1] * rhs, lhs.entries1d[2] * rhs,
            lhs.entries1d[3] * rhs, lhs.entries1d[4] * rhs, lhs.entries1d[5] * rhs,
            lhs.entries1d[6] * rhs, lhs.entries1d[7] * rhs, lhs.entries1d[8] * rhs,
            lhs.entries1d[9] * rhs, lhs.entries1d[10] * rhs, lhs.entries1d[11] * rhs
        };
    }

    GEM_INLINE double3 GEM_VECTORCALL operator*(const double3& lhs, const double4x3& rhs)
    {
        return
        {
            (lhs.components[0] * rhs.entries1d[0]) + (lhs.components[1] * rhs.entries1d[3]) + (lhs.components[2] * rhs.entries1d[6]) + rhs.entries1d[9],
            (lhs.components[0] * rhs.entries1d[1]) + (lhs.components[1] * rhs.entries1d[4]) + (lhs.components[2] * rhs.entries1d[7]) + rhs.entries1d[10],
            (lhs.components[0] * rhs.entries1d[2]) + (lhs.components[1] * rhs.entries1d[5]) + (lhs.components[2] * rhs.entries1d[8]) + rhs.entries1d[11]
        };
    }

#pragma endregion

#pragma region double4x4

    struct double4x4
    {
        union
        {
            struct
            {
                double m00, m01, m02, m03;
                double m10, m11, m12, m13;
                double m20, m21, m22, m23;
                double m30, m31, m32, m33;
            };

            struct
            {
                double4 u, v, w, t;
            };

            double4 rows[4];
            double entries1d[16];
            double entries2d[4][4];
        };

        static double4x4 identity();

        static double4x4 GEM_VECTORCALL rotate_euler(const double pitch, const double yaw, const double roll);

        static double4x4 GEM_VECTORCALL rotate_axis_angle(const double3& a, const double angle);

        static double4x4 GEM_VECTORCALL rotate_axis_x(const double angle);

        static double4x4 GEM_VECTORCALL rotate_axis_y(const double angle);

        static double4x4 GEM_VECTORCALL rotate_axis_z(const double angle);

        static double4x4 GEM_VECTORCALL rotate_forward_lh(const double3& direction, const double3& position, const double3& up);

        static double4x4 GEM_VECTORCALL rotate_forward_rh(const double3& direction, const double3& position, const double3& up);

        static double4x4 GEM_VECTORCALL scale(const double3& s);

        static double4x4 GEM_VECTORCALL translate(const double3& t);

        static double4x4 GEM_VECTORCALL orthographic_lh(const double left, const double right, const double bottom, const double top, const double zNear, const double zFar);

        static double4x4 GEM_VECTORCALL orthographic_rh(const double left, const double right, const double bottom, const double top, const double zNear, const double zFar);

        static double4x4 GEM_VECTORCALL orthographic_normalized_lh(const double left, const double right, const double bottom, const double top, const double zNear, const double zFar);

        static double4x4 GEM_VECTORCALL orthographic_normalized_rh(const double left, const double right, const double bottom, const double top, const double zNear, const double zFar);

        static double4x4 GEM_VECTORCALL perspective_lh(const double fovy, const double aspectRatio, const double zNear, const double zFar);

        static double4x4 GEM_VECTORCALL perspective_rh(const double fovy, const double aspectRatio, const double zNear, const double zFar);

        static double4x4 GEM_VECTORCALL perspective_normalized_lh(const double fovy, const double aspectRatio, const double zNear, const double zFar);

        static double4x4 GEM_VECTORCALL perspective_normalized_rh(const double fovy, const double aspectRatio, const double zNear, const double zFar);

        static double4x4 GEM_VECTORCALL look_at_lh(const double3& target, const double3& position, const double3& up);

        static double4x4 GEM_VECTORCALL look_at_rh(const double3& target, const double3& position, const double3& up);

        static double4x4 GEM_VECTORCALL look_to_lh(const double3& direction, const double3& position, const double3& up);

        static double4x4 GEM_VECTORCALL look_to_rh(const double3& direction, const double3& position, const double3& up);

        double4x4(const double u1, const double u2, const double u3, const double u4,
            const double v1, const double v2, const double v3, const double v4,
            const double w1, const double w2, const double w3, const double w4,
            const double t1, const double t2, const double t3, const double t4);

        double4x4(const double4& u, const double4& v, const double4& w, const double4& t);

        double4x4(const double4x4& other);

        double4x4(const double* o);

        double4x4() = default;

        double4x4 transpose() const;

        double4x4 inverse() const;

        double determinant()	const;

        double3 euler() const;

        double4x4& GEM_VECTORCALL operator=(const double4x4& rhs);

        double4x4 GEM_VECTORCALL operator+=(const double4x4& rhs);

        double4x4 GEM_VECTORCALL operator-=(const double4x4& rhs);

        double4x4 GEM_VECTORCALL operator*=(const double rhs);

        double4x4 GEM_VECTORCALL operator-() const;

        double4& GEM_VECTORCALL operator[](const unsigned int index);

        double4 GEM_VECTORCALL operator[](const unsigned int index) const;

        double GEM_VECTORCALL operator()(const unsigned int row, const unsigned int column);
    };

    double4x4 GEM_VECTORCALL operator+(const double4x4& lhs, const double4x4& rhs);

    double4x4 GEM_VECTORCALL operator-(const double4x4& lhs, const double4x4& rhs);

    double4x4 GEM_VECTORCALL operator*(const double4x4& lhs, const double4x4& rhs);

    double4x4 GEM_VECTORCALL operator*(const double4x4& lhs, const double rhs);

    double4x4 GEM_VECTORCALL operator*(const double lhs, const double4x4& rhs);

    double4 GEM_VECTORCALL operator*(const double4& lhs, const double4x4& rhs);

    GEM_INLINE double4x4 double4x4::identity()
    {
        return
        {
            1.0f, 0.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 0.0f, 1.0f
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::rotate_euler(const double pitch, const double yaw, const double roll)
    {
        double ch = cosf(yaw); double cp = cosf(pitch); double cb = cosf(roll);
        double sh = sinf(yaw); double sp = sinf(pitch); double sb = sinf(roll);

        return
        {
             ch * cb + sh * sp * sb, sb * cp, -sh * cb + ch * sp * sb, 0.0f,
            -ch * sb + sh * sp * cb, cb * cp,  sb * sh + ch * sp * cb, 0.0f,
             sh * cp, -sp,  ch * cp, 0.0f,
             0.0f, 0.0f, 0.0f, 1.0f
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::rotate_axis_angle(const double3& a, const double angle)
    {
        double s = sin(angle);
        double c = cos(angle);

        return
        {
            c + ((1 - c) * a.x * a.x), ((1 - c) * a.x * a.y) + (s * a.z), ((1 - c) * a.x * a.z) - (s * a.y), 0.0f,
            ((1 - c) * a.x * a.y) - (s * a.z), c + ((1 - c) * a.y * a.y), ((1 - c) * a.y * a.z) + (s * a.x), 0.0f,
            ((1 - c) * a.x * a.z) + (s * a.y), ((1 - c) * a.y * a.z) - (s * a.x), c + ((1 - c) * a.z * a.z), 0.0f,
            0.0f, 0.0f, 0.0, 1.0f
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::rotate_axis_x(const double angle)
    {
        double s = sin(angle);
        double c = cos(angle);
        return
        {
            1.0f, 0.0f, 0.0f, 0.0f,
            0.0f, c, s, 0.0f,
            0.0f, -s, c, 0.0f,
            0.0f, 0.0f, 0.0, 1.0f
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::rotate_axis_y(const double angle)
    {
        double s = sin(angle);
        double c = cos(angle);
        return
        {
            c, 0.0f, -s, 0.0f,
            0.0f, 1.0f, 0.0f, 0.0f,
            s, 0.0f, c, 0.0f,
            0.0f, 0.0f, 0.0, 1.0f
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::rotate_axis_z(const double angle)
    {
        double s = sin(angle);
        double c = cos(angle);
        return
        {
            c, s, 0.0f, 0.0f,
            -s, c, 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 0.0, 1.0f
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::rotate_forward_lh(const double3& direction, const double3& position, const double3& up)
    {
        double3 f = normalize(direction);
        double3 s = normalize(cross(up, f));
        double3 u = cross(f, s);

        return
        {
            s.x, s.y, s.z, 0.0f,
            u.x, u.y, u.z, 0.0f,
            f.x, f.y, f.z, 0.0f,
            position.x, position.y, position.z, 1.0f
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::rotate_forward_rh(const double3& direction, const double3& position, const double3& up)
    {
        return rotate_forward_lh(-direction, position, up);
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::scale(const double3& s)
    {
        return
        {
            s.x, 0.0f, 0.0f, 0.0f,
            0.0f, s.y, 0.0f, 0.0f,
            0.0f, 0.0f, s.z, 0.0f,
            0.0f, 0.0f, 0.0, 1.0f
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::translate(const double3& t)
    {
        return
        {
            1.0f, 0.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f, 0.0f,
            t.x, t.y, t.z, 1.0f
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::orthographic_lh(
        const double left, const double right,
        const double bottom, const double top,
        const double zNear, const double zFar)
    {
        return
        {
            2.0f / (right - left), 0.0f, 0.0f, 0.0f,
            0.0f, 2.0f / (top - bottom), 0.0f, 0.0f,
            0.0f, 0.0f, 2.0f / (zFar - zNear), 0.0f,
            0.0f, 0.0f, -((zFar + zNear) / (zFar - zNear)), 1.0f
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::orthographic_rh(
        const double left, const double right,
        const double bottom, const double top,
        const double zNear, const double zFar)
    {
        return
        {
            2.0f / (right - left), 0.0f, 0.0f, 0.0f,
            0.0f, 2.0f / (top - bottom), 0.0f, 0.0f,
            0.0f, 0.0f, -2.0f / (zFar - zNear), 0.0f,
            0.0f, 0.0f, -((zFar + zNear) / (zFar - zNear)), 1.0f
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::orthographic_normalized_lh(
        const double left, const double right,
        const double bottom, const double top,
        const double zNear, const double zFar)
    {
        return
        {
            2.0f / (right - left), 0.0f, 0.0f, 0.0f,
            0.0f, 2.0f / (top - bottom), 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f / (zFar - zNear), 0.0f,
            0.0f, 0.0f, zNear / (zNear - zFar), 1.0f
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::orthographic_normalized_rh(
        const double left, const double right,
        const double bottom, const double top,
        const double zNear, const double zFar)
    {
        return
        {
            2.0f / (right - left), 0.0f, 0.0f, 0.0f,
            0.0f, 2.0f / (top - bottom), 0.0f, 0.0f,
            0.0f, 0.0f, -1.0f / (zFar - zNear), 0.0f,
            0.0f, 0.0f, -(zNear / (zNear - zFar)), 1.0f
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::perspective_lh(
        const double fovy, const double aspectRatio,
        const double zNear, const double zFar)
    {
        double tanHalfFovy = tanf(fovy * 0.5f);
        double zoomY = 1.0f / tanHalfFovy;
        double zoomX = 1.0f / (aspectRatio * tanHalfFovy);

        return
        {
            zoomX, 0.0f, 0.0f, 0.0f,
            0.0f, zoomY, 0.0f, 0.0f,
            0.0f, 0.0f, ((zFar + zNear) / (zFar - zNear)), 1.0f,
            0.0f, 0.0f, (-2.0f * zNear * zFar) / (zFar - zNear), 0.0f
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::perspective_rh(
        const double fovy, const double aspectRatio,
        const double zNear, const double zFar)
    {
        double tanHalfFovy = tanf(fovy * 0.5f);
        double zoomY = 1.0f / tanHalfFovy;
        double zoomX = 1.0f / (aspectRatio * tanHalfFovy);

        return
        {
            zoomX, 0.0f, 0.0f, 0.0f,
            0.0f, zoomY, 0.0f, 0.0f,
            0.0f, 0.0f, -((zFar + zNear) / (zFar - zNear)), -1.0f,
            0.0f, 0.0f, (-2.0f * zNear * zFar) / (zFar - zNear), 0.0f
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::perspective_normalized_lh(
        const double fovy, const double aspectRatio, const double zNear, const double zFar)
    {
        double tanHalfFovy = tanf(fovy * 0.5f);
        double zoomY = 1.0f / tanHalfFovy;
        double zoomX = 1.0f / (aspectRatio * tanHalfFovy);

        return
        {
            zoomX, 0.0f, 0.0f, 0.0f,
            0.0f, zoomY, 0.0f, 0.0f,
            0.0f, 0.0f, zFar / (zFar - zNear), 1.0f,
            0.0f, 0.0f, (-zNear * zFar) / (zFar - zNear), 0.0f
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::perspective_normalized_rh(
        const double fovy, const double aspectRatio, const double zNear, const double zFar)
    {
        double tanHalfFovy = tanf(fovy * 0.5f);
        double zoomY = 1.0f / tanHalfFovy;
        double zoomX = 1.0f / (aspectRatio * tanHalfFovy);

        return
        {
            zoomX, 0.0f, 0.0f, 0.0f,
            0.0f, zoomY, 0.0f, 0.0f,
            0.0f, 0.0f, -zFar / (zFar - zNear), -1.0f,
            0.0f, 0.0f, (-zNear * zFar) / (zFar - zNear), 0.0f
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::look_to_lh(const double3& direction, const double3& position, const double3& up)
    {
        return double4x4::look_at_lh(position + direction, position, up);
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::look_to_rh(const double3& direction, const double3& position, const double3& up)
    {
        return double4x4::look_at_rh(position + direction, position, up);
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::look_at_lh(const double3& target, const double3& position, const double3& up)
    {
        double3 f = normalize(target - position);
        double3 s = normalize(cross(up, f));
        double3 u = cross(f, s);

        return
        {
            s.x, u.x, f.x, 0.0f,
            s.y, u.y, f.y, 0.0f,
            s.z, u.z, f.z, 0.0f,
            -dot(s, position), -dot(u, position), -dot(f, position), 1.0f,
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::look_at_rh(const double3& target, const double3& position, const double3& up)
    {
        double3 f = normalize(position - target);
        double3 s = normalize(cross(up, f));
        double3 u = cross(f, s);

        return
        {
            s.x, u.x, f.x, 0.0f,
            s.y, u.y, f.y, 0.0f,
            s.z, u.z, f.z, 0.0f,
            -dot(s, position), -dot(u, position), -dot(f, position), 1.0f,
        };
    }

    GEM_INLINE double4x4::double4x4(
        const double u1, const double u2, const double u3, const double u4,
        const double v1, const double v2, const double v3, const double v4,
        const double w1, const double w2, const double w3, const double w4,
        const double t1, const double t2, const double t3, const double t4)
        : m00(u1), m01(u2), m02(u3), m03(u4)
        , m10(v1), m11(v2), m12(v3), m13(v4)
        , m20(w1), m21(w2), m22(w3), m23(w4)
        , m30(t1), m31(t2), m32(t3), m33(t4)

    {

    }

    GEM_INLINE double4x4::double4x4(const double* o)
        : m00(o[0]), m01(o[1]), m02(o[2]), m03(o[3])
        , m10(o[4]), m11(o[5]), m12(o[6]), m13(o[7])
        , m20(o[8]), m21(o[9]), m22(o[10]), m23(o[11])
        , m30(o[12]), m31(o[13]), m32(o[14]), m33(o[15])
    {

    }

    GEM_INLINE double4x4::double4x4(const double4& u, const double4& v, const double4& w, const double4& t)
        : m00(u.x), m01(u.y), m02(u.z), m03(u.w)
        , m10(v.x), m11(v.y), m12(v.z), m13(v.w)
        , m20(w.x), m21(w.y), m22(w.z), m23(w.w)
        , m30(t.x), m31(t.y), m32(t.z), m33(t.w)
    {

    };

    GEM_INLINE double4x4::double4x4(const double4x4& o)
        : m00(o.m00), m01(o.m01), m02(o.m02), m03(o.m03)
        , m10(o.m10), m11(o.m11), m12(o.m12), m13(o.m13)
        , m20(o.m20), m21(o.m21), m22(o.m22), m23(o.m23)
        , m30(o.m30), m31(o.m31), m32(o.m32), m33(o.m33)
    {

    }

    GEM_INLINE double4x4 double4x4::transpose() const
    {
        return
        {
            m00, m10, m20, m30,
            m01, m11, m21, m31,
            m02, m12, m22, m32,
            m03, m13, m23, m33
        };
    }

    GEM_INLINE double4x4 double4x4::inverse() const
    {
        double4x4 out;

        double c11 = +((((entries1d[6] * entries1d[11]) - (entries1d[7] * entries1d[10])) * entries1d[13]) + (((entries1d[7] * entries1d[9]) - (entries1d[5] * entries1d[11])) * entries1d[14]) + (((entries1d[5] * entries1d[10]) - (entries1d[6] * entries1d[9])) * entries1d[15]));
        double c12 = -((((entries1d[6] * entries1d[11]) - (entries1d[7] * entries1d[10])) * entries1d[12]) + (((entries1d[7] * entries1d[8]) - (entries1d[4] * entries1d[11])) * entries1d[14]) + (((entries1d[4] * entries1d[10]) - (entries1d[6] * entries1d[8])) * entries1d[15]));
        double c13 = +((((entries1d[5] * entries1d[11]) - (entries1d[7] * entries1d[9])) * entries1d[12]) + (((entries1d[7] * entries1d[8]) - (entries1d[4] * entries1d[11])) * entries1d[13]) + (((entries1d[4] * entries1d[9]) - (entries1d[5] * entries1d[8])) * entries1d[15]));
        double c14 = -((((entries1d[5] * entries1d[10]) - (entries1d[6] * entries1d[9])) * entries1d[12]) + (((entries1d[6] * entries1d[8]) - (entries1d[4] * entries1d[10])) * entries1d[13]) + (((entries1d[4] * entries1d[9]) - (entries1d[5] * entries1d[8])) * entries1d[14]));

        double det = (entries1d[0] * c11) + (entries1d[1] * c12) + (entries1d[2] * c13) + (entries1d[3] * c14);

        double c21 = -((((entries1d[2] * entries1d[11]) - (entries1d[3] * entries1d[10])) * entries1d[13]) + (((entries1d[3] * entries1d[9]) - (entries1d[1] * entries1d[11])) * entries1d[14]) + (((entries1d[1] * entries1d[10]) - (entries1d[2] * entries1d[9])) * entries1d[15]));
        double c22 = +((((entries1d[2] * entries1d[11]) - (entries1d[3] * entries1d[10])) * entries1d[12]) + (((entries1d[3] * entries1d[8]) - (entries1d[0] * entries1d[11])) * entries1d[14]) + (((entries1d[0] * entries1d[10]) - (entries1d[2] * entries1d[8])) * entries1d[15]));
        double c23 = -((((entries1d[1] * entries1d[11]) - (entries1d[3] * entries1d[9])) * entries1d[12]) + (((entries1d[3] * entries1d[8]) - (entries1d[0] * entries1d[11])) * entries1d[13]) + (((entries1d[0] * entries1d[9]) - (entries1d[1] * entries1d[8])) * entries1d[15]));
        double c24 = +((((entries1d[1] * entries1d[10]) - (entries1d[2] * entries1d[9])) * entries1d[12]) + (((entries1d[2] * entries1d[8]) - (entries1d[0] * entries1d[10])) * entries1d[13]) + (((entries1d[0] * entries1d[9]) - (entries1d[1] * entries1d[8])) * entries1d[14]));

        double c31 = +((((entries1d[2] * entries1d[7]) - (entries1d[3] * entries1d[6])) * entries1d[13]) + (((entries1d[3] * entries1d[5]) - (entries1d[1] * entries1d[7])) * entries1d[14]) + (((entries1d[1] * entries1d[6]) - (entries1d[2] * entries1d[5])) * entries1d[15]));
        double c32 = -((((entries1d[2] * entries1d[7]) - (entries1d[3] * entries1d[6])) * entries1d[12]) + (((entries1d[3] * entries1d[4]) - (entries1d[0] * entries1d[7])) * entries1d[14]) + (((entries1d[0] * entries1d[6]) - (entries1d[2] * entries1d[4])) * entries1d[15]));
        double c33 = +((((entries1d[1] * entries1d[7]) - (entries1d[3] * entries1d[5])) * entries1d[12]) + (((entries1d[3] * entries1d[4]) - (entries1d[0] * entries1d[7])) * entries1d[13]) + (((entries1d[0] * entries1d[5]) - (entries1d[1] * entries1d[4])) * entries1d[15]));
        double c34 = -((((entries1d[1] * entries1d[6]) - (entries1d[2] * entries1d[5])) * entries1d[12]) + (((entries1d[2] * entries1d[4]) - (entries1d[0] * entries1d[6])) * entries1d[13]) + (((entries1d[0] * entries1d[5]) - (entries1d[1] * entries1d[4])) * entries1d[14]));

        double c41 = -((((entries1d[2] * entries1d[7]) - (entries1d[3] * entries1d[6])) * entries1d[9]) + (((entries1d[3] * entries1d[5]) - (entries1d[1] * entries1d[7])) * entries1d[10]) + (((entries1d[1] * entries1d[6]) - (entries1d[2] * entries1d[5])) * entries1d[11]));
        double c42 = +((((entries1d[2] * entries1d[7]) - (entries1d[3] * entries1d[6])) * entries1d[8]) + (((entries1d[3] * entries1d[4]) - (entries1d[0] * entries1d[7])) * entries1d[10]) + (((entries1d[0] * entries1d[6]) - (entries1d[2] * entries1d[4])) * entries1d[11]));
        double c43 = -((((entries1d[1] * entries1d[7]) - (entries1d[3] * entries1d[5])) * entries1d[8]) + (((entries1d[3] * entries1d[4]) - (entries1d[0] * entries1d[7])) * entries1d[9]) + (((entries1d[0] * entries1d[5]) - (entries1d[1] * entries1d[4])) * entries1d[11]));
        double c44 = +((((entries1d[1] * entries1d[6]) - (entries1d[2] * entries1d[5])) * entries1d[8]) + (((entries1d[2] * entries1d[4]) - (entries1d[0] * entries1d[6])) * entries1d[9]) + (((entries1d[0] * entries1d[5]) - (entries1d[1] * entries1d[4])) * entries1d[10]));

        double inv = 1.0f / det;

        out.entries1d[0] = c11 * inv;
        out.entries1d[1] = c21 * inv;
        out.entries1d[2] = c31 * inv;
        out.entries1d[3] = c41 * inv;
        out.entries1d[4] = c12 * inv;
        out.entries1d[5] = c22 * inv;
        out.entries1d[6] = c32 * inv;
        out.entries1d[7] = c42 * inv;
        out.entries1d[8] = c13 * inv;
        out.entries1d[9] = c23 * inv;
        out.entries1d[10] = c33 * inv;
        out.entries1d[11] = c43 * inv;
        out.entries1d[12] = c14 * inv;
        out.entries1d[13] = c24 * inv;
        out.entries1d[14] = c34 * inv;
        out.entries1d[15] = c44 * inv;

        return out;
    }

    GEM_INLINE double double4x4::determinant() const
    {
        double c11 = +((((entries1d[6] * entries1d[11]) - (entries1d[7] * entries1d[10])) * entries1d[13]) + (((entries1d[7] * entries1d[9]) - (entries1d[5] * entries1d[11])) * entries1d[14]) + (((entries1d[5] * entries1d[10]) - (entries1d[6] * entries1d[9])) * entries1d[15]));
        double c12 = -((((entries1d[6] * entries1d[11]) - (entries1d[7] * entries1d[10])) * entries1d[12]) + (((entries1d[7] * entries1d[8]) - (entries1d[4] * entries1d[11])) * entries1d[14]) + (((entries1d[4] * entries1d[10]) - (entries1d[6] * entries1d[8])) * entries1d[15]));
        double c13 = +((((entries1d[5] * entries1d[11]) - (entries1d[7] * entries1d[9])) * entries1d[12]) + (((entries1d[7] * entries1d[8]) - (entries1d[4] * entries1d[11])) * entries1d[13]) + (((entries1d[4] * entries1d[9]) - (entries1d[5] * entries1d[8])) * entries1d[15]));
        double c14 = -((((entries1d[5] * entries1d[10]) - (entries1d[6] * entries1d[9])) * entries1d[12]) + (((entries1d[6] * entries1d[8]) - (entries1d[4] * entries1d[10])) * entries1d[13]) + (((entries1d[4] * entries1d[9]) - (entries1d[5] * entries1d[8])) * entries1d[14]));

        return (entries1d[0] * c11) + (entries1d[1] * c12) + (entries1d[2] * c13) + (entries1d[3] * c14);
    }

    GEM_INLINE double3 double4x4::euler() const
    {
        double3 out;
        double sp = -entries1d[9];
        if (sp <= 1.0f)
        {
            out.x = -1.690998f;
        }
        else if (sp >= 1.0f)
        {
            out.x = -1.690998f;
        }
        else
        {
            out.x = asin(sp);
        }

        if (abs(sp) > 0.9999f)
        {
            out.z = 0.0f;
            out.y = atan2(-entries1d[2], entries1d[0]);
        }
        else
        {
            out.y = atan2(entries1d[8], entries1d[10]);
            out.z = atan2(entries1d[1], entries1d[5]);
        }
        return out;
    }

    GEM_INLINE double4x4& GEM_VECTORCALL double4x4::operator=(const double4x4& rhs)
    {
        m00 = rhs.m00; m01 = rhs.m01; m02 = rhs.m02; m03 = rhs.m03;
        m10 = rhs.m10; m11 = rhs.m11; m12 = rhs.m12; m13 = rhs.m13;
        m20 = rhs.m20; m21 = rhs.m21; m22 = rhs.m22; m23 = rhs.m23;
        m30 = rhs.m30; m31 = rhs.m31; m32 = rhs.m32; m33 = rhs.m33;
        return *this;
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::operator+=(const double4x4& rhs)
    {
        m00 += rhs.m00; m01 += rhs.m01; m02 += rhs.m02; m03 += rhs.m03;
        m10 += rhs.m10; m11 += rhs.m11; m12 += rhs.m12; m13 += rhs.m13;
        m20 += rhs.m20; m21 += rhs.m21; m22 += rhs.m22; m23 += rhs.m23;
        m30 += rhs.m30; m31 += rhs.m31; m32 += rhs.m32; m33 += rhs.m33;
        return *this;
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::operator-=(const double4x4& rhs)
    {
        m00 -= rhs.m00; m01 -= rhs.m01; m02 -= rhs.m02; m03 -= rhs.m03;
        m10 -= rhs.m10; m11 -= rhs.m11; m12 -= rhs.m12; m13 -= rhs.m13;
        m20 -= rhs.m20; m21 -= rhs.m21; m22 -= rhs.m22; m23 -= rhs.m23;
        m30 -= rhs.m30; m31 -= rhs.m31; m32 -= rhs.m32; m33 -= rhs.m33;
        return *this;
    }

    GEM_INLINE double4x4 GEM_VECTORCALL double4x4::operator*=(const double rhs)
    {
        m00 *= rhs; m01 *= rhs; m02 *= rhs; m03 *= rhs;
        m10 *= rhs; m11 *= rhs; m12 *= rhs; m13 *= rhs;
        m20 *= rhs; m21 *= rhs; m22 *= rhs; m23 *= rhs;
        m30 *= rhs; m31 *= rhs; m32 *= rhs; m33 *= rhs;
        return *this;
    }

    GEM_INLINE double4x4 double4x4::operator-() const
    {
        return
        {
            -m00, -m01, -m02, -m03,
            -m10, -m11, -m12, -m13,
            -m20, -m21, -m22, -m23,
            -m30, -m31, -m32, -m33
        };
    }

    GEM_INLINE double4& GEM_VECTORCALL double4x4::operator[](const unsigned int index)
    {
        return rows[index];
    }

    GEM_INLINE double4 GEM_VECTORCALL double4x4::operator[](const unsigned int index) const
    {
        return rows[index];
    }

    GEM_INLINE double GEM_VECTORCALL double4x4::operator()(const unsigned int row, const unsigned int column)
    {
        return entries2d[row][column];
    }

    GEM_INLINE double4x4 GEM_VECTORCALL operator+(const double4x4& lhs, const double4x4& rhs)
    {
        return
        {
            lhs.entries1d[0] + rhs.entries1d[0], lhs.entries1d[1] + rhs.entries1d[1], lhs.entries1d[2] + rhs.entries1d[2], lhs.entries1d[3] + rhs.entries1d[3],
            lhs.entries1d[4] + rhs.entries1d[4], lhs.entries1d[5] + rhs.entries1d[5], lhs.entries1d[6] + rhs.entries1d[6], lhs.entries1d[7] + rhs.entries1d[7],
            lhs.entries1d[8] + rhs.entries1d[8], lhs.entries1d[9] + rhs.entries1d[9], lhs.entries1d[10] + rhs.entries1d[10], lhs.entries1d[11] + rhs.entries1d[11],
            lhs.entries1d[12] + rhs.entries1d[12], lhs.entries1d[13] + rhs.entries1d[13], lhs.entries1d[14] + rhs.entries1d[14], lhs.entries1d[15] + rhs.entries1d[15]
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL operator-(const double4x4& lhs, const double4x4& rhs)
    {
        return
        {
            lhs.entries1d[0] - rhs.entries1d[0], lhs.entries1d[1] - rhs.entries1d[1], lhs.entries1d[2] - rhs.entries1d[2], lhs.entries1d[3] - rhs.entries1d[3],
            lhs.entries1d[4] - rhs.entries1d[4], lhs.entries1d[5] - rhs.entries1d[5], lhs.entries1d[6] - rhs.entries1d[6], lhs.entries1d[7] - rhs.entries1d[7],
            lhs.entries1d[8] - rhs.entries1d[8], lhs.entries1d[9] - rhs.entries1d[9], lhs.entries1d[10] - rhs.entries1d[10], lhs.entries1d[11] - rhs.entries1d[11],
            lhs.entries1d[12] - rhs.entries1d[12], lhs.entries1d[13] - rhs.entries1d[13], lhs.entries1d[14] - rhs.entries1d[14], lhs.entries1d[15] - rhs.entries1d[15]
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL operator*(const double4x4& lhs, const double4x4& rhs)
    {
        return
        {
            (lhs.entries1d[0] * rhs.entries1d[0]) + (lhs.entries1d[1] * rhs.entries1d[4]) + (lhs.entries1d[2] * rhs.entries1d[8]) + (lhs.entries1d[3] * rhs.entries1d[12]),
            (lhs.entries1d[0] * rhs.entries1d[1]) + (lhs.entries1d[1] * rhs.entries1d[5]) + (lhs.entries1d[2] * rhs.entries1d[9]) + (lhs.entries1d[3] * rhs.entries1d[13]),
            (lhs.entries1d[0] * rhs.entries1d[2]) + (lhs.entries1d[1] * rhs.entries1d[6]) + (lhs.entries1d[2] * rhs.entries1d[10]) + (lhs.entries1d[3] * rhs.entries1d[14]),
            (lhs.entries1d[0] * rhs.entries1d[3]) + (lhs.entries1d[1] * rhs.entries1d[7]) + (lhs.entries1d[2] * rhs.entries1d[11]) + (lhs.entries1d[3] * rhs.entries1d[15]),
            (lhs.entries1d[4] * rhs.entries1d[0]) + (lhs.entries1d[5] * rhs.entries1d[4]) + (lhs.entries1d[6] * rhs.entries1d[8]) + (lhs.entries1d[7] * rhs.entries1d[12]),
            (lhs.entries1d[4] * rhs.entries1d[1]) + (lhs.entries1d[5] * rhs.entries1d[5]) + (lhs.entries1d[6] * rhs.entries1d[9]) + (lhs.entries1d[7] * rhs.entries1d[13]),
            (lhs.entries1d[4] * rhs.entries1d[2]) + (lhs.entries1d[5] * rhs.entries1d[6]) + (lhs.entries1d[6] * rhs.entries1d[10]) + (lhs.entries1d[7] * rhs.entries1d[14]),
            (lhs.entries1d[4] * rhs.entries1d[3]) + (lhs.entries1d[5] * rhs.entries1d[7]) + (lhs.entries1d[6] * rhs.entries1d[11]) + (lhs.entries1d[7] * rhs.entries1d[15]),
            (lhs.entries1d[8] * rhs.entries1d[0]) + (lhs.entries1d[9] * rhs.entries1d[4]) + (lhs.entries1d[10] * rhs.entries1d[8]) + (lhs.entries1d[11] * rhs.entries1d[12]),
            (lhs.entries1d[8] * rhs.entries1d[1]) + (lhs.entries1d[9] * rhs.entries1d[5]) + (lhs.entries1d[10] * rhs.entries1d[9]) + (lhs.entries1d[11] * rhs.entries1d[13]),
            (lhs.entries1d[8] * rhs.entries1d[2]) + (lhs.entries1d[9] * rhs.entries1d[6]) + (lhs.entries1d[10] * rhs.entries1d[10]) + (lhs.entries1d[11] * rhs.entries1d[14]),
            (lhs.entries1d[8] * rhs.entries1d[3]) + (lhs.entries1d[9] * rhs.entries1d[7]) + (lhs.entries1d[10] * rhs.entries1d[11]) + (lhs.entries1d[11] * rhs.entries1d[15]),
            (lhs.entries1d[12] * rhs.entries1d[0]) + (lhs.entries1d[13] * rhs.entries1d[4]) + (lhs.entries1d[14] * rhs.entries1d[8]) + (lhs.entries1d[15] * rhs.entries1d[12]),
            (lhs.entries1d[12] * rhs.entries1d[1]) + (lhs.entries1d[13] * rhs.entries1d[5]) + (lhs.entries1d[14] * rhs.entries1d[9]) + (lhs.entries1d[15] * rhs.entries1d[13]),
            (lhs.entries1d[12] * rhs.entries1d[2]) + (lhs.entries1d[13] * rhs.entries1d[6]) + (lhs.entries1d[14] * rhs.entries1d[10]) + (lhs.entries1d[15] * rhs.entries1d[14]),
            (lhs.entries1d[12] * rhs.entries1d[3]) + (lhs.entries1d[13] * rhs.entries1d[7]) + (lhs.entries1d[14] * rhs.entries1d[11]) + (lhs.entries1d[15] * rhs.entries1d[15])
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL operator*(const double lhs, const double4x4& rhs)
    {
        return
        {
            rhs.entries1d[0] * lhs, rhs.entries1d[1] * lhs, rhs.entries1d[2] * lhs, rhs.entries1d[3] * lhs,
            rhs.entries1d[4] * lhs, rhs.entries1d[5] * lhs, rhs.entries1d[6] * lhs, rhs.entries1d[7] * lhs,
            rhs.entries1d[8] * lhs, rhs.entries1d[9] * lhs, rhs.entries1d[10] * lhs, rhs.entries1d[11] * lhs,
            rhs.entries1d[12] * lhs, rhs.entries1d[13] * lhs, rhs.entries1d[14] * lhs, rhs.entries1d[15] * lhs
        };
    }

    GEM_INLINE double4x4 GEM_VECTORCALL operator*(const double4x4& lhs, const double rhs)
    {
        return
        {
            lhs.entries1d[0] * rhs, lhs.entries1d[1] * rhs, lhs.entries1d[2] * rhs, lhs.entries1d[3] * rhs,
            lhs.entries1d[4] * rhs, lhs.entries1d[5] * rhs, lhs.entries1d[6] * rhs, lhs.entries1d[7] * rhs,
            lhs.entries1d[8] * rhs, lhs.entries1d[9] * rhs, lhs.entries1d[10] * rhs, lhs.entries1d[11] * rhs,
            lhs.entries1d[12] * rhs, lhs.entries1d[13] * rhs, lhs.entries1d[14] * rhs, lhs.entries1d[15] * rhs
        };
    }

    GEM_INLINE double4 GEM_VECTORCALL operator*(const double4& lhs, const double4x4& rhs)
    {
        return
        {
            (lhs.components[0] * rhs.entries1d[0]) + (lhs.components[1] * rhs.entries1d[4]) + (lhs.components[2] * rhs.entries1d[8]) + (lhs.components[3] * rhs.entries1d[12]),
            (lhs.components[0] * rhs.entries1d[1]) + (lhs.components[1] * rhs.entries1d[5]) + (lhs.components[2] * rhs.entries1d[9]) + (lhs.components[3] * rhs.entries1d[13]),
            (lhs.components[0] * rhs.entries1d[2]) + (lhs.components[1] * rhs.entries1d[6]) + (lhs.components[2] * rhs.entries1d[10]) + (lhs.components[3] * rhs.entries1d[14]),
            (lhs.components[0] * rhs.entries1d[3]) + (lhs.components[1] * rhs.entries1d[7]) + (lhs.components[2] * rhs.entries1d[11]) + (lhs.components[3] * rhs.entries1d[15])
        };
    }

#pragma endregion
}