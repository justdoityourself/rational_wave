/* Copyright (C) 2020 D8DATAWORKS - All Rights Reserved */

#pragma once

#include "rw/rational_wave.hpp"

using namespace rational_wave;

TEST_CASE("RW Addition", "[rw::]")
{
	RationalWave<int16_t> r1(Repeat(5), 5,8,3,1);
	RationalWave<int16_t> r2(Repeat(4),9,3,7,2,1);
	RationalWave<int16_t> sum(Repeat(1), 14, 11, 10, 3, 6, 17, 6, 8, 7, 9, 12, 4, 12, 10, 4, 10, 8, 15, 5, 2);

	r2 += r1;

	CHECK((r2 == sum));
}

TEST_CASE("Frequency Isolation", "[rw::]")
{
	RationalWave<int16_t> r1(Repeat(5), 5, 8, 3, 1);
	RationalWave<int16_t> r2(Repeat(4), 9, 3, 7, 2, 1);
	RationalWave<int16_t> sum(Repeat(1), 14, 11, 10, 3, 6, 17, 6, 8, 7, 9, 12, 4, 12, 10, 4, 10, 8, 15, 5, 2);

	auto freq4 = sum.ApplyTransformation(DerivativeIdentity<int16_t>(5 % 4, 4));
	auto _freq5 = sum.ApplyTransformation(DerivativeIdentity<int16_t>(4,5));
	auto freq5 = _freq5.AntiDerivative(sum[0] / 2);

	CHECK((freq4 == r1));
	CHECK((freq5 == r2));
}

TEST_CASE("Legacy", "[rw::]")
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
