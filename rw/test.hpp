/* Copyright (C) 2020 D8DATAWORKS - All Rights Reserved */

#pragma once

#include "rw/rational_wave.hpp"

using namespace rational_wave;


TEST_CASE("Basic", "[rw::]")
{

	/*
		Add two waves of different frequency, and separate them again... validating the results;
	*/

	DerivativeIdentity<int16_t> t(7, 13);

	t.Print();

	RationalWave<int16_t> r1(19, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10);
	RationalWave<int16_t> r2(11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18);

	r2 += r1;

	auto dr = r2.DerivativeFrequency(19);
	auto br = dr.IsBalancedAndRegular(19);

	auto mr = r2.DerivativeFrequency(11);
	br = mr.IsBalancedAndRegular(11);

	auto mrprime = mr.ApplyTransformation(DerivativeIdentity<int16_t>(19 % 11, 11));

	auto tr = dr.ApplyTransformation(DerivativeIdentity<int16_t>(11, 19));
	//dr.Invert();

	auto final = tr.AntiDerivative(r2[0] / 2);
	br = final.IsBalancedAndRegular(19);

	auto final2 = mrprime.AntiDerivative(r2[0] / 2);
	br = final2.IsBalancedAndRegular(11);

	r2.Console();
	final.Console();
	final2.Console();
}
