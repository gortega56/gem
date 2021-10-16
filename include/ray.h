#pragma once
#include "plane.h"
#include "range.h"
#include "transform.h"

namespace gem
{
    struct ray3f
    {
        float3 p;
        float3 v;
    
        bool intersects_ray(const ray3f& ray, float3* p_point, const float tolerance = 0.001f) const;

        bool intersects_plane(const plane4f& plane, float3* p_point, const float tolerance = 0.001f) const;

        bool intersects_range(const range3f& range, float3* p_point, const float tolerance = 0.001f) const;
    };

    bool ray3f::intersects_ray(const ray3f& ray, float3* p_point, const float tolerance /*= 0.001f*/) const
    {
        float3 v0 = v;
        float3 v1 = ray.v;
        float3 p0 = p;
        float3 p1 = ray.p;
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

    bool ray3f::intersects_plane(const plane4f& plane, float3* p_point, const float tolerance /*= 0.001f*/) const
    {
        float ndotv = dot(plane.n, v);
        if (ndotv < tolerance)
            return false;

        if (p_point)
            *p_point = p - v * (dot(plane.n, p) / ndotv);

        return true;
    }

    bool ray3f::intersects_range(const range3f& range, float3* p_point, const float tolerance /*= 0.001f*/) const
    {
        https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
        float3 iv =
        {
            1.0f / v.x,
            1.0f / v.y,
            1.0f / v.z
        };

        int sign[3] =
        {
            iv.x < 0.0f,
            iv.y < 0.0f,
            iv.z < 0.0f
        };

        float3 bounds[2] =
        {
            range.min,
            range.max
        };

        float tmin, tmax, tymin, tymax, tzmin, tzmax;
        tmin = (bounds[sign[0]].x - p.x) * iv.x;
        tmax = (bounds[1 - sign[0]].x - p.x) * iv.x;
        tymin = (bounds[sign[1]].y - p.y) * iv.y;
        tymax = (bounds[1 - sign[1]].y - p.y) * iv.y;

        if ((tmin > tymax) || (tymin > tmax))
            return false;
        if (tymin > tmin)
            tmin = tymin;
        if (tymax < tmax)
            tmax = tymax;

        tzmin = (bounds[sign[2]].z - p.z) * iv.z;
        tzmax = (bounds[1 - sign[2]].z - p.z) * iv.z;

        if ((tmin > tzmax) || (tzmin > tmax))
            return false;
        if (tzmin > tmin)
            tmin = tzmin;
        if (tzmax < tmax)
            tmax = tzmax;

        if (p_point)
            *p_point = p + tmin * v;

        return true;
    }
}