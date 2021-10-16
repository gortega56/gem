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

    }

    bool ray3f::intersects_plane(const plane4f& plane, float3* p_point, const float tolerance /*= 0.001f*/) const
    {

    }

    bool ray3f::intersects_range(const range3f& range, float3* p_point, const float tolerance /*= 0.001f*/) const
    {

    }
}