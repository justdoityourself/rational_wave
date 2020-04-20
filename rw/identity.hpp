/* Copyright (C) 2020 D8DATAWORKS - All Rights Reserved */

#pragma once

#include <vector>
#include <iostream>

namespace rational_wave
{
    /*
        A DerivativeIdentify takes two frequencies and generates a matrix solution which isolates each waveform.
    */

    template < typename T > class DerivativeIdentity
    {
        struct N
        {
            N(T* _p) :p(_p) {}

            T& operator[](T i) const
            {
                return *(p + i);
            }

            T* p;
        };
    public:
        DerivativeIdentity(T _x, T _y)
            : x(_x)
            , y(_y)
        {
            if (x > y)
            {
                T t = x;
                x = y;
                y = t;
            }

            std::vector<T> tmp;
            tmp.resize(y * y);
            T* p = tmp.data();

            auto abs = [](T v)-> T { if (v > 0) return v; return v * -1; };

            T d = y - x;
            T positive = x;
            T negative = 1;

            T pitr = positive;
            T nitr = negative;

            T itr = 2;
            p[0] = positive;
            p[1] = negative;
            while (d != abs(pitr - nitr))
            {
                pitr = (pitr + x) % y;

                nitr = nitr - x;
                if (nitr < 0)
                    nitr += y;

                p[itr + 1] = nitr;
                p[itr] = pitr;

                itr += 2;
            }

            sl = itr;

            m.resize(y * sl);
            T* mp = m.data();

            for (size_t i = 0; i < y; i++)
                for (size_t j = 0; j < sl; j++)
                    *mp++ = (p[j] + i) % y;
        }

        N operator[](T i) const
        {
            return N(((T*)m.data()) + i * sl);
        }

        void Print()
        {
            for (size_t i = 0; i < y; i++)
            {
                for (size_t j = 0; j < x; j++)
                    std::cout << (int64_t)(*this)[i][j] << " ";

                std::cout << std::endl;
            }
        }

        size_t height() const
        {
            return y;
        }

        size_t length() const
        {
            return sl;
        }

    private:
        size_t sl;
        size_t x;
        size_t y;
        std::vector<T> m;
    };
}