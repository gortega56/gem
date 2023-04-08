#pragma once
#include "transform.h"

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

    struct line3f
    {
        float3 v, m;

        float3 GEM_VECTORCALL point_closest_to_line(const line3f& line, const float tolerance = 0.001f) const;

        float GEM_VECTORCALL distance_to_line(const line3f& line, const float tolerance = 0.001f) const;

        float3 GEM_VECTORCALL point_closest_to_point(const float3& point) const;

        float GEM_VECTORCALL distance_to_point(const float3& point) const;

        bool GEM_VECTORCALL intersects_line(const line3f& line, float3* p_point, const float tolerance = 0.001f) const;

        line3f& GEM_VECTORCALL transform(const transform3f& transform);

        line3f& GEM_VECTORCALL transform(const transform3f& transform);
    };

    GEM_INLINE float3 GEM_VECTORCALL line3f::point_closest_to_line(const line3f& line, const float tolerance /*= 0.001f*/) const
    {
        float3 v0 = v;
        float3 v1 = line.v;
        float3 p0 = cross(v0, m) * (1.0f / length_squared(v0));
        float3 p1 = cross(v1, line.m) * (1.0f / length_squared(v1));
        float3 dp = p1 - p0;
        float v00 = length_squared(v0);
        float v11 = length_squared(v1);
        float v01 = dot(v0, v1);
        float det = (v01 * v01) - (v00 * v11);
        if (det < tolerance)
            return sqrtf(length_squared(cross(dp, v1)) / v01);

        det = 1.0f / det;
        float dpv0 = dot(dp, v0);
        float dpv1 = dot(dp, v1);
        float t = ((v01 * dpv1) - (v11 * dpv0)) * det;
        return p0 + v0 * t;
    }

    GEM_INLINE float GEM_VECTORCALL line3f::distance_to_line(const line3f& line, const float tolerance /*= 0.001f*/) const
    {
        return (dot(v, line.m) + dot(line.v, m)) / length(cross(v, line.v));
    }

    GEM_INLINE float3 GEM_VECTORCALL line3f::point_closest_to_point(const float3& point) const
    {
        float3 closest_to_origin = cross(v, m) * (1.0f / length_squared(v));
        return project(closest_to_origin - point, v);
    }

    GEM_INLINE float GEM_VECTORCALL line3f::distance_to_point(const float3& point) const
    {
        return length(cross(v, point) + m) / length(v);
    }

    GEM_INLINE bool GEM_VECTORCALL line3f::intersects_line(const line3f& line, float3* p_point, const float tolerance /*= 0.001f*/) const
    {
        float3 v0 = v;
        float3 v1 = line.v;
        float3 p0 = cross(v0, m) * (1.0f / length_squared(v0));
        float3 p1 = cross(v1, line.m) * (1.0f / length_squared(v1));
        float3 dp = p1 - p0;
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

        if (p_point)
            *p_point = p0 + v0 * t;

        return true;
    }

    GEM_INLINE line3f& GEM_VECTORCALL line3f::transform(const transform3f& transform)
    {
        float4x3 h = transform.matrix4x3().inverse().transpose();
        v = transform.transform_vector(v);
        m = m * h;
        return *this;
    }

    GEM_INLINE line3f& GEM_VECTORCALL line3f::transform(const transform3f& transform)
    {
        float4x3 h = transform.matrix4x3().inverse().transpose();
        v = transform.transform_vector(v);
        m = m * h;
        return *this;
    }

    struct line_segment2f
    {
        float2 a, b;

        float2 GEM_VECTORCALL point_closest_to_point(const float2& point) const;
        
        float GEM_VECTORCALL distance_to_point(const float2& point) const;

        bool GEM_VECTORCALL intersects_line(const line_segment2f& line, float2* p_point, const float tolerance = 0.001f) const;
    };

    GEM_INLINE float2 GEM_VECTORCALL line_segment2f::point_closest_to_point(const float2& point) const
    {
        return project(point - a, b - a);
    }

    GEM_INLINE float GEM_VECTORCALL line_segment2f::distance_to_point(const float2& point) const
    {
        float d = (point.x * a.y) - (point.y * a.x);
        float l = length_squared(b - a);
        return sqrtf((d * d) / l);
    }

    GEM_INLINE bool GEM_VECTORCALL line_segment2f::intersects_line(const line_segment2f& line, float2* p_point, const float tolerance /*= 0.001f*/) const
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

    struct line_segment3f
    {
        float3 a, b;

        float3 GEM_VECTORCALL point_closest_to_line(const line_segment3f& line, const float tolerance = 0.001f) const;
        
        float GEM_VECTORCALL distance_to_line(const line_segment3f& line, const float tolerance = 0.001f) const;

        float3 GEM_VECTORCALL point_closest_to_point(const float3& point) const;
        
        float GEM_VECTORCALL distance_to_point(const float3& point) const;

        bool GEM_VECTORCALL intersects_line(const line_segment3f& line, float3* p_point, const float tolerance = 0.001f) const;

        line_segment3f& GEM_VECTORCALL transform(const transform3f& transform);

        line_segment3f& GEM_VECTORCALL transform(const transform3f& transform);
    };

    GEM_INLINE float3 GEM_VECTORCALL line_segment3f::point_closest_to_line(const line_segment3f& line, const float tolerance /*= 0.001f*/) const
    {
        float3 v0 = b - a;
        float3 v1 = line.b - line.a;
        float3 dp = line.a - a;
        float v00 = length_squared(v0);
        float v11 = length_squared(v1);
        float v01 = dot(v0, v1);
        float det = (v01 * v01) - (v00 * v11);
        if (det < tolerance)
            return sqrtf(length_squared(cross(dp, v1)) / v01);

        det = 1.0f / det;
        float dpv0 = dot(dp, v0);
        float dpv1 = dot(dp, v1);
        float t = ((v01 * dpv1) - (v11 * dpv0)) * det;
        if (t < 0) t = 0;
        if (t > 1) t = 1;
        return line.a + v0 * t;
    }

    GEM_INLINE float GEM_VECTORCALL line_segment3f::distance_to_line(const line_segment3f& line, const float tolerance /*= 0.001f*/) const
    {
        float3 v0 = b - a;
        float3 v1 = line.b - line.a;
        float3 dp = line.a - a;
        float v00 = length_squared(v0);
        float v11 = length_squared(v1);
        float v01 = dot(v0, v1);
        float det = (v01 * v01) - (v00 * v11);
        if (det < tolerance)
            return sqrtf(length_squared(cross(dp, v1)) / v01);

        det = 1.0f / det;
        float dpv0 = dot(dp, v0);
        float dpv1 = dot(dp, v1);
        float t0 = ((v01 * dpv1) - (v11 * dpv0)) * det;
        float t1 = ((v00 * dpv1) - (v01 * dpv0)) * det;
        if (t0 < 0) t0 = 0;
        if (t0 > 1) t0 = 1;
        if (t1 < 0) t1 = 0;
        if (t1 > 1) t1 = 1;
        return length(dp + (v1 * t1) - (v0 * t0));
    }

    GEM_INLINE float3 GEM_VECTORCALL line_segment3f::point_closest_to_point(const float3& point) const
    {
        return project(point - a, b - a);
    }

    GEM_INLINE float GEM_VECTORCALL line_segment3f::distance_to_point(const float3& point) const
    {
        float d = length_squared(cross(point, a));
        float l = length_squared(b - a);
        return sqrtf(d / l);
    }

    GEM_INLINE bool GEM_VECTORCALL line_segment3f::intersects_line(const line_segment3f& line, float3* p_point, const float tolerance /*= 0.001f*/) const
    {
        float3 v0 = b - a;
        float3 v1 = line.b - line.a;
        float3 dp = line.a - a;
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

    GEM_INLINE line_segment3f& GEM_VECTORCALL line_segment3f::transform(const transform3f& transform)
    {
        a = transform.transform_point(a);
        b = transform.transform_point(b);
        return *this;
    }

    GEM_INLINE line_segment3f& GEM_VECTORCALL line_segment3f::transform(const transform3f& transform)
    {
        a = transform.transform_point(a);
        b = transform.transform_point(b);
        return *this;
    }
}