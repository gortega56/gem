#include <gtest/gtest.h>

#include <functional>
#include <random>

#include "../include/vector.h"
#include "../include/matrix.h"
#include "../include/quaternion.h"
#include "../include/transform.h"

#include "matrix.h"
#include "vector.h"
#include "quaternion.h"

#define GEM_ASSERT 1

namespace gem
{
	namespace unit_tests
	{
		struct float1
		{
			float1(float components)
				: x(components)
			{

			}

			float1() = default;

			union
			{
				struct
				{
					float x;
				};

				float components[1];
			};
		};

		typedef float1 control_float1;
		typedef lina::vector<float, 2> control_float2;
		typedef lina::vector<float, 3> control_float3;
		typedef lina::vector<float, 4> control_float4;
		
		typedef lina::matrix<float, 2, 2> control_float2x2;
		typedef lina::matrix<float, 3, 3> control_float3x3;
		typedef lina::matrix<float, 4, 3> control_float4x3;
		typedef lina::matrix<float, 4, 4> control_float4x4;

		typedef lina::quaternion<float> control_quatf;

		typedef float1 vary_float1;
		typedef float2 vary_float2;
		typedef float3 vary_float3;
		typedef float4 vary_float4;

		typedef float2x2 vary_float2x2;
		typedef float3x3 vary_float3x3;
		typedef float4x3 vary_float4x3;
		typedef float4x4 vary_float4x4;

		typedef quatf vary_quatf;

		static std::random_device rd;
		static std::mt19937 gen(rd());

		template<typename return_type, typename operand_type>
		struct unary_op
		{
			typedef return_type(function)(operand_type&);
			return_type result;
			operand_type a;
			function* p_fn;
			unary_op(function* fn)
				: result()
				, a()
				, p_fn(fn)
			{

			}

			void operator()()
			{
				result = p_fn(a);
			}
		};

		template<typename return_type, typename operand_type0, typename operand_type1>
		struct binary_op
		{
			typedef return_type(function)(operand_type0&, operand_type1&);
			return_type result;
			operand_type0 a;
			operand_type1 b;
			function* p_fn;
			binary_op(function* fn)
				: result()
				, a()
				, b()
				, p_fn(fn)
			{

			}

			void operator()()
			{
				result = p_fn(a, b);
			}
		};


		template<typename T, unsigned int num_components>
		void fill(std::mt19937& gen, T* p_expected, T* p_actual)
		{
			std::uniform_real_distribution<> dis(-10.0f, 10.0f);
			for (unsigned int i = 0; i < num_components; ++i)
				p_expected[i] = p_actual[i] = static_cast<T>(dis(gen));
		}

		template<typename T, unsigned int num_components>
		void expect_near(T* p_expected, T* p_actual, float tolerance)
		{
			for (unsigned int i = 0; i < num_components; ++i)
				EXPECT_NEAR(p_expected[i], p_actual[i], tolerance);
		}

		template<typename component_type, typename expected_type, typename actual_type, unsigned int num_components_return_type, unsigned int num_components_operand_type0, unsigned int num_components_operand_type1>
		void test_binary_op(std::mt19937& gen, expected_type& expected, actual_type& actual, float tolerance)
		{
			fill<component_type, num_components_operand_type0>(gen, reinterpret_cast<component_type*>(&expected.a), reinterpret_cast<component_type*>(&actual.a));
			fill<component_type, num_components_operand_type1>(gen, reinterpret_cast<component_type*>(&expected.b), reinterpret_cast<component_type*>(&actual.b));
			expected();
			actual();
			expect_near<component_type, num_components_return_type>(reinterpret_cast<component_type*>(&expected.result), reinterpret_cast<component_type*>(&actual.result), tolerance);
		}

		template<typename component_type, typename expected_type, typename actual_type, unsigned int num_components_return_type, unsigned int num_components_operand_type>
		void test_unary_op(std::mt19937& gen, expected_type& expected, actual_type& actual, float tolerance)
		{
			fill<component_type, num_components_operand_type>(gen, reinterpret_cast<component_type*>(&expected.a), reinterpret_cast<component_type*>(&actual.a));
			expected();
			actual();
			expect_near<component_type, num_components_return_type>(reinterpret_cast<component_type*>(&expected.result), reinterpret_cast<component_type*>(&actual.result), tolerance);
		}

#define RUN_UNARY_OPERATION_CASE0(underlying_type, num_components_return_type, num_components_operand_type, tolerance, operation)\
			typedef control_ ## underlying_type ## num_components_return_type  control_return_type;	\
			typedef control_ ## underlying_type ## num_components_operand_type control_operand_type;\
			typedef vary_ ## underlying_type ## num_components_return_type vary_return_type;\
			typedef vary_ ## underlying_type ## num_components_operand_type vary_operand_type;\
			typedef unary_op<control_return_type, control_operand_type> control_op;\
			typedef unary_op<vary_return_type, vary_operand_type> vary_op;\
			control_op expected([](control_operand_type& a) { return (operation); });\
			vary_op actual([](vary_operand_type& a) { return (operation); });\
			test_unary_op<underlying_type, control_op, vary_op, num_components_return_type, num_components_operand_type>(gen, expected, actual, tolerance)

#define RUN_UNARY_OPERATION_CASE1(underlying_type, return_type, operand_type, num_components_return_type, num_components_operand_type, tolerance, operation)\
			typedef control_ ## return_type control_return_type;\
			typedef control_ ## operand_type control_operand_type;\
			typedef vary_ ## return_type vary_return_type;\
			typedef vary_ ## operand_type vary_operand_type;\
			typedef unary_op<control_return_type, control_operand_type> control_op;\
			typedef unary_op<vary_return_type, vary_operand_type> vary_op;\
			control_op expected([](control_operand_type& a) { return (operation); });\
			vary_op actual([](vary_operand_type& a) { return (operation); });\
			test_unary_op<underlying_type, control_op, vary_op, num_components_return_type, num_components_operand_type>(gen, expected, actual, tolerance)

#define RUN_BINARY_OPERATION_CASE0(underlying_type, num_components_return_type, num_components_operand_type0, num_components_operand_type1, tolerance, operation)\
			typedef control_ ## underlying_type ## num_components_return_type	control_return_type;\
			typedef control_ ## underlying_type ## num_components_operand_type0 control_operand_type0;\
			typedef control_ ## underlying_type ## num_components_operand_type1 control_operand_type1;\
			typedef vary_ ## underlying_type ## num_components_return_type vary_return_type;\
			typedef vary_ ## underlying_type ## num_components_operand_type0 vary_operand_type0;\
			typedef vary_ ## underlying_type ## num_components_operand_type0 vary_operand_type1;\
			typedef binary_op<control_return_type, control_operand_type0, control_operand_type1> control_op;\
			typedef binary_op<vary_return_type, vary_operand_type0, vary_operand_type1> vary_op;\
			control_op expected([](control_operand_type0& a, control_operand_type1& b) { return (operation); });\
			vary_op actual([](vary_operand_type0& a, vary_operand_type1& b) { return (operation); });\
			test_binary_op<underlying_type, control_op, vary_op, num_components_return_type, num_components_operand_type0, num_components_operand_type1>(gen, expected, actual, tolerance)

#define RUN_BINARY_OPERATION_CASE1(underlying_type, return_type, operand_type0, operand_type1, num_components_return_type, num_components_operand_type0, num_components_operand_type1, tolerance, operation)\
			typedef control_ ## return_type control_return_type;\
			typedef control_ ## operand_type0 control_operand_type0;\
			typedef control_ ## operand_type1 control_operand_type1;\
			typedef vary_ ## return_type vary_return_type;\
			typedef vary_ ## operand_type0 vary_operand_type0;\
			typedef vary_ ## operand_type1 vary_operand_type1;\
			typedef binary_op<control_return_type, control_operand_type0, control_operand_type1> control_op;\
			typedef binary_op<vary_return_type, vary_operand_type0, vary_operand_type1> vary_op;\
			control_op expected([](control_operand_type0& a, control_operand_type1& b) { return (operation); });\
			vary_op actual([](vary_operand_type0& a, vary_operand_type1& b) { return (operation); });\
			test_binary_op<underlying_type, control_op, vary_op, num_components_return_type, num_components_operand_type0, num_components_operand_type1>(gen, expected, actual, tolerance)

#define EXPECT_NEAR_ARRAY(underlying_type, num_components, expected, actual, tolerance) expect_near<underlying_type, num_components>(reinterpret_cast<underlying_type*>(&expected), reinterpret_cast<underlying_type*>(&actual), tolerance)
		TEST(float2, addition)
		{
			RUN_BINARY_OPERATION_CASE0(float, 2, 2, 2, 0.001f, a + b);
		}

		TEST(float2, subtraction)
		{
			RUN_BINARY_OPERATION_CASE0(float, 2, 2, 2, 0.001f, a - b);
		}

		TEST(float2, multiplication)
		{
			RUN_BINARY_OPERATION_CASE0(float, 2, 2, 2, 0.001f, a * b);
		}

		TEST(float2, dot)
		{
			RUN_BINARY_OPERATION_CASE0(float, 1, 2, 2, 0.001f, float1(dot(a,b)));
		}

		TEST(float2, projection)
		{
			RUN_BINARY_OPERATION_CASE0(float, 2, 2, 2, 0.001f, project(a, b) + reject(a,b));
		}

		TEST(float2, length_squared)
		{
			RUN_UNARY_OPERATION_CASE0(float, 1, 2, 0.001f, float1(length_squared(a)));
		}

		TEST(float2, length)
		{
			RUN_UNARY_OPERATION_CASE0(float, 1, 2, 0.001f, float1(length(a)));
		}

		TEST(float3, addition)
		{
			RUN_BINARY_OPERATION_CASE0(float, 3, 3, 3, 0.001f, a + b);
		}

		TEST(float3, subtraction)
		{
			RUN_BINARY_OPERATION_CASE0(float, 3, 3, 3, 0.001f, a - b);
		}

		TEST(float3, multiplication)
		{
			RUN_BINARY_OPERATION_CASE0(float, 3, 3, 3, 0.001f, a * b);
		}

		TEST(float3, dot)
		{
			RUN_BINARY_OPERATION_CASE0(float, 1, 3, 3, 0.001f, float1(dot(a, b)));
		}

		TEST(float3, length_squared)
		{
			RUN_UNARY_OPERATION_CASE0(float, 1, 3, 0.001f, float1(length_squared(a)));
		}

		TEST(float3, length)
		{
			RUN_UNARY_OPERATION_CASE0(float, 1, 3, 0.001f, float1(length(a)));
		}

		TEST(float4, addition)
		{
			RUN_BINARY_OPERATION_CASE0(float, 4, 4, 4, 0.001f, a + b);
		}

		TEST(float4, subtraction)
		{
			RUN_BINARY_OPERATION_CASE0(float, 4, 4, 4, 0.001f, a - b);
		}

		TEST(float4, multiplication)
		{
			RUN_BINARY_OPERATION_CASE0(float, 4, 4, 4, 0.001f, a * b);
		}

		TEST(float4, dot)
		{
			RUN_BINARY_OPERATION_CASE0(float, 1, 4, 4, 0.001f, float1(dot(a, b)));
		}

		TEST(float4, length_squared)
		{
			RUN_UNARY_OPERATION_CASE0(float, 1, 4, 0.001f, float1(length_squared(a)));
		}

		TEST(float4, length)
		{
			RUN_UNARY_OPERATION_CASE0(float, 1, 4, 0.001f, float1(length(a)));
		}

		TEST(float2x2, addition)
		{
			RUN_BINARY_OPERATION_CASE1(float, float2x2, float2x2, float2x2, 4, 4, 4, 0.001f, a + b);
		}

		TEST(float2x2, subtraction)
		{
			RUN_BINARY_OPERATION_CASE1(float, float2x2, float2x2, float2x2, 4, 4, 4, 0.001f, a - b);
		}

		TEST(float2x2, multiplication)
		{
			RUN_BINARY_OPERATION_CASE1(float, float2x2, float2x2, float2x2, 4, 4, 4, 0.001f, a * b);
		}

		TEST(float2x2, transformation)
		{
			RUN_BINARY_OPERATION_CASE1(float, float2, float2, float2x2, 2, 4, 4, 0.001f, a * b);
		}

		TEST(float2x2, determinant)
		{
			RUN_UNARY_OPERATION_CASE1(float, float1, float2x2, 1, 4, 0.001f, float1(a.determinant()));
		}

		TEST(float2x2, inverse)
		{
			RUN_UNARY_OPERATION_CASE1(float, float2x2, float2x2, 1, 4, 0.001f, a * a.inverse());
		}

		TEST(float2x2, transpose)
		{
			RUN_UNARY_OPERATION_CASE1(float, float2x2, float2x2, 1, 4, 0.001f, a.transpose());
		}

		TEST(float3x3, addition)
		{
			RUN_BINARY_OPERATION_CASE1(float, float3x3, float3x3, float3x3, 9, 9, 9, 0.001f, a + b);
		}

		TEST(float3x3, subtraction)
		{
			RUN_BINARY_OPERATION_CASE1(float, float3x3, float3x3, float3x3, 9, 9, 9, 0.001f, a - b);
		}

		TEST(float3x3, multiplication)
		{
			RUN_BINARY_OPERATION_CASE1(float, float3x3, float3x3, float3x3, 9, 9, 9, 0.001f, a * b);
		}

		TEST(float3x3, transformation)
		{
			RUN_BINARY_OPERATION_CASE1(float, float3, float3, float3x3, 3, 9, 9, 0.001f, a * b);
		}

		TEST(float3x3, determinant)
		{
			RUN_UNARY_OPERATION_CASE1(float, float1, float3x3, 1, 9, 0.001f, float1(a.determinant()));
		}

		TEST(float3x3, inverse)
		{
			RUN_UNARY_OPERATION_CASE1(float, float3x3, float3x3, 1, 9, 0.001f, a * a.inverse());
		}

		TEST(float3x3, transpose)
		{
			RUN_UNARY_OPERATION_CASE1(float, float3x3, float3x3, 1, 9, 0.001f, a.transpose());
		}

		TEST(float4x4, addition)
		{
			RUN_BINARY_OPERATION_CASE1(float, float4x4, float4x4, float4x4, 16, 16, 16, 0.001f, a + b);
		}

		TEST(float4x4, subtraction)
		{
			RUN_BINARY_OPERATION_CASE1(float, float4x4, float4x4, float4x4, 16, 16, 16, 0.001f, a - b);
		}

		TEST(float4x4, multiplication)
		{
			RUN_BINARY_OPERATION_CASE1(float, float4x4, float4x4, float4x4, 16, 16, 16, 0.001f, a * b);
		}

		TEST(float4x4, transformation)
		{
			RUN_BINARY_OPERATION_CASE1(float, float4, float4, float4x4, 3, 16, 16, 0.001f, a * b);
		}

		TEST(float4x4, determinant)
		{
			RUN_UNARY_OPERATION_CASE1(float, float1, float4x4, 1, 16, 0.001f, float1(a.determinant()));
		}

		TEST(float4x4, inverse)
		{
			RUN_UNARY_OPERATION_CASE1(float, float4x4, float4x4, 1, 16, 0.001f, a * a.inverse());
		}

		TEST(float4x4, transpose)
		{
			RUN_UNARY_OPERATION_CASE1(float, float4x4, float4x4, 1, 16, 0.001f, a.transpose());
		}

		TEST(quatf, multiplication)
		{
			RUN_BINARY_OPERATION_CASE1(float, quatf, quatf, quatf, 4, 4, 4, 0.001f, a * b);
		}

		TEST(quatf, transformation)
		{
			RUN_BINARY_OPERATION_CASE1(float, float3, float3, quatf, 3, 3, 4, 0.001f, b * a);
		}

		TEST(transform3f, concatenation)
		{
			float3 euler = float3(0.261799f, 0.261799f, 0.261799f);
			float3 scale = float3(0.3f, 0.4f, 1.2f);
			float3 xlate = float3(-10.0f, 59.0f, -2.0f);

			float4x3 expected;
			{
				float4x3 q = float4x3::rotate_euler(euler.x, euler.y, euler.z);
				float4x3 s = float4x3::scale(scale);
				float4x3 t = float4x3::translate(xlate);
				expected = (q * s * t);
			}

			float4x3 actual;
			{
				transform3f q = transform3f::set(quatf::rotate_euler(euler.x, euler.y, euler.z), float3(0, 0, 0), float3(1, 1, 1));
				transform3f s = transform3f::set(quatf::identity(), float3(0, 0, 0), scale);
				transform3f t = transform3f::set(quatf::identity(), xlate, float3(1, 1, 1));
				transform3f xform = transform3f::identity();
				xform.concatenate(q);
				xform.concatenate(s);
				xform.concatenate(t);
				actual = xform.matrix4x3();
			}

			EXPECT_NEAR_ARRAY(float, 12, expected, actual, 0.001f);
		}

		TEST(transform3f, transformation)
		{
			float3 euler = float3(0.261799f, 0.261799f, 0.261799f);
			float3 scale = float3(0.3f, 0.4f, 1.2f);
			float3 xlate = float3(-10.0f, 59.0f, -2.0f);
			float3 pos = float3(-73.0f, 0.5f, 14.8f);

			float3 expected;
			{
				float4x3 q = float4x3::rotate_euler(euler.x, euler.y, euler.z);
				float4x3 s = float4x3::scale(scale);
				float4x3 t = float4x3::translate(xlate);
				expected = pos * (q * s * t);
			}

			float3 actual;
			{
				transform3f xform = transform3f::set(quatf::rotate_euler(euler.x, euler.y, euler.z), xlate, scale);
				actual = xform.transform_point(pos);
			}

			EXPECT_NEAR_ARRAY(float, 3, expected, actual, 0.001f);
		}

		TEST(transform1f, concatenation)
		{
			float  scale = 1.2f;
			float3 euler = float3(0.261799f, 0.261799f, 0.261799f);
			float3 xlate = float3(-10.0f, 59.0f, -2.0f);

			float4x3 expected;
			{
				float4x3 q = float4x3::rotate_euler(euler.x, euler.y, euler.z);
				float4x3 s = float4x3::scale(scale);
				float4x3 t = float4x3::translate(xlate);
				expected = (q * s * t);
			}

			float4x3 actual;
			{
				transform1f q = transform1f::set(quatf::rotate_euler(euler.x, euler.y, euler.z), float3(0, 0, 0), 1);
				transform1f s = transform1f::set(quatf::identity(), float3(0, 0, 0), scale);
				transform1f t = transform1f::set(quatf::identity(), xlate, 1);
				transform1f xform = transform1f::identity();
				xform.concatenate(q);
				xform.concatenate(s);
				xform.concatenate(t);
				actual = xform.matrix4x3();
			}

			EXPECT_NEAR_ARRAY(float, 12, expected, actual, 0.001f);
		}

		TEST(transform1f, transformation)
		{
			float  scale = 1.2f;
			float3 euler = float3(0.261799f, 0.261799f, 0.261799f);
			float3 xlate = float3(-10.0f, 59.0f, -2.0f);
			float3 pos = float3(-73.0f, 0.5f, 14.8f);

			float3 expected;
			{
				float4x3 q = float4x3::rotate_euler(euler.x, euler.y, euler.z);
				float4x3 s = float4x3::scale(scale);
				float4x3 t = float4x3::translate(xlate);
				expected = pos * (q * s * t);
			}

			float3 actual;
			{
				transform1f xform = transform1f::set(quatf::rotate_euler(euler.x, euler.y, euler.z), xlate, scale);
				actual = xform.transform_point(pos);
			}

			EXPECT_NEAR_ARRAY(float, 3, expected, actual, 0.001f);
		}
	}
}