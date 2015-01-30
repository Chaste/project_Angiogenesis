/*

Copyright (c) 2005-2015, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#ifndef TESTGEOMETRYTRANSFORM_HPP_
#define TESTGEOMETRYTRANSFORM_HPP_

#include <cxxtest/TestSuite.h>
#include "SmartPointers.hpp"
#include "FakePetscSetup.hpp"
#include "ChastePoint.hpp"
#include "GeometryTransform.hpp"

class TestGeometryTransform : public CxxTest::TestSuite
{
public:

	void TestTranslating() throw(Exception)
	{
		// Make some points
		ChastePoint<2> point1(2.0, 3.0);
		ChastePoint<2> point2(-4.0, 3.5);
		std::vector<ChastePoint<2> > points_2d;
		points_2d.push_back(point1);
		points_2d.push_back(point2);

		ChastePoint<3> point3(2.0, 3.0, 4.0);
		ChastePoint<3> point4(-4.0, 3.5, 8.8);
		std::vector<ChastePoint<3> > points_3d;
		points_3d.push_back(point3);
		points_3d.push_back(point4);

		// Move them
		std::vector<double> translation_vector_2d;
		translation_vector_2d.push_back(3.5);
		translation_vector_2d.push_back(5.6);

		std::vector<double> translation_vector_3d = translation_vector_2d;
		translation_vector_3d.push_back(-12.8);

		GeometryTransform<2> transform_2d;
		GeometryTransform<3> transform_3d;

		std::vector<ChastePoint<2> > new_points_2d = transform_2d.Translate(points_2d, translation_vector_2d);
		std::vector<ChastePoint<3> > new_points_3d = transform_3d.Translate(points_3d, translation_vector_3d);

		// Check some of the resulting co-ordinates
		TS_ASSERT_DELTA(new_points_3d[0][0], 5.5, 1.e-6);
		TS_ASSERT_DELTA(new_points_2d[1][0], -0.5, 1.e-6);
		TS_ASSERT_DELTA(new_points_3d[1][2], -4.0, 1.e-6);
	}
};

#endif /*TESTGEOMETRYTRANSFORMHPP_*/
