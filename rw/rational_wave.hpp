/* Copyright (C) 2020 D8DATAWORKS - All Rights Reserved */

#pragma once

#include <iostream>

#include "identity.hpp"

namespace rational_wave
{
    /*
        A Rational Wave wraps an sample buffer and provides methods to comprehend the data.
    */

    class Repeat
    {
    public:
        Repeat(size_t _r) : repetition(_r) {}

        size_t operator() () const
        {
            return repetition;
        }

    private:
        size_t repetition;
    };

    template <typename T > class RationalWave
    {
    public:
        RationalWave(){}

        template <typename ... F> RationalWave(Repeat rep, F ... args)
        {
            Generate(rep,args...);
        }

        template <uint64_t O> void _Generate(){}

        template <uint64_t O, typename Q, typename ... F> void _Generate(Q c, F ... args)
        {
            (*this)[O] = (T)c;
            _Generate<O+1>(args...);
        }

        template <typename ... F> void Generate(Repeat rep, F ... args)
        {
            constexpr size_t unit =  sizeof...(args);
            set_size(unit*rep());

            _Generate<0>(args...);

            T gap = unit;
            for(size_t i = 1; i < rep(); i++)
                std::copy(m.begin(),m.begin()+gap,m.data() + gap * i);
        }

        template < typename F > void PairIterator(const RationalWave & r, F f)
        {
            if(r.m.size() != m.size())
                throw 1;

            T * p1 = data();
            const T * p2 = r.data_ro();

            size_t l = m.size();
            while(l--)
                f(*p1++,*p2++);		
        }

        void operator -=(const RationalWave & r)
        {
            PairIterator(r,[](T & l, const T & r) { l -= r; });
        }

        void operator +=(const RationalWave & r)
        {
            PairIterator(r,[](T & l, const T & r) { l += r; });	
        }

        bool operator==(const RationalWave& r)
        {
            if (r.samples() != samples())
                return false;

            return std::equal(m.begin(), m.end(), r.m.begin());
        }

        T & operator[](size_t i){ return m[i]; }

        size_t samples() const { return m.size(); }
        const T* data_ro() const { return m.data(); }
        T * data() { return m.data(); }

        void set_size(size_t c)
        {
            m.resize(c);
        }

        RationalWave DerivativeFrequency(T frequency)
        {
            RationalWave rw;

            size_t l = m.size() / frequency;

            if(m.size() % frequency)
                throw std::runtime_error("This function requires exact multiple");

            size_t sz = samples();
            rw.set_size(sz);

            for(size_t i = 0; i<sz;i++)
                rw[(i + l) % sz] = (*this)[(i + l) % sz] - (*this)[i];

            return rw;
        }

        void Multiply(T m)
        {
            auto r = data();
            auto s = samples();

            for(T i = 0;i<s;i++)
                r[i] *= m;
        }

        void Invert()
        {
            Multiply(-1);
        }

        RationalWave ApplyTransformation(const DerivativeIdentity<T> & d)
        {
            RationalWave rw;

            size_t sz = samples();
            if(sz % d.height())
                throw std::runtime_error("Wrong Identity");

            size_t r = sz / d.height();

            rw.set_size(sz);

            for(size_t i = 0; i < d.height(); i++)
            {
                for(size_t j = 0; j < d.length(); j++)
                {
                    if(j==0)
                        rw[i] = (*this)[d[i][j]];
                    else
                        rw[i] += (*this)[d[i][j]];
                }
            }

            for(size_t i = 1; i < r; i++ )
                std::copy(rw.m.begin(), rw.m.begin() + d.height(), rw.m.data() + i * d.height());

            return rw;
        }

        RationalWave AntiDerivative(T amplitude)
        {
            RationalWave rw;

            size_t sz = samples();
            rw.set_size(sz);

            rw[0] = amplitude;

            for(size_t i = 1; i < sz; i++)
            {
                amplitude += (*this)[i-1];
                rw[i] = amplitude;
            }

            return rw;
        }

        std::pair<bool,bool> IsBalancedAndRegular(T frequency)
        {
            bool balance = true;
            bool regular = true;

            size_t l = m.size() / frequency;
            if(m.size() % frequency)
                throw std::runtime_error("This function requires exact multiple");

            size_t sz = samples();
            T sum = 0;

            for(size_t i = 0; i < sz; i++)
            {
                if((*this)[i] != (*this)[(i + frequency)%sz])
                    regular = false;

                if(i % frequency == 0)
                {
                    if(sum)
                    {
                        balance = false;
                        sum = 0;
                    }
                }
                sum += (*this)[i];
            }

            if(sum)
            {
                balance = false;
                sum = 0;
            }

            return std::make_pair(balance,regular);
        }

        RationalWave BalanceWave()
        {
            RationalWave rw;

            return rw;
        }

        void Console()
        {
            auto r = data();
            
            size_t l = samples();

            for(size_t i = 0; i < l;i++)
                std::cout << (int64_t)r[i] << " ";
            
            std::cout << std::endl << std::endl;
        }

    private:
        std::vector<T> m;
    };
}