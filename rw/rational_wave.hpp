#include <vector>
#include <iostream>

namespace rational_wave
{
	#define self_t (*this)

    /*
        TODO Relocate || Refactor this object.
    */

    class ManagedMemory : public std::vector<uint8_t>
    {
    public:
        ManagedMemory(){}
        /*ManagedMemory(const Memory & m )
        {
            resize(m.size());
            std::memcpy(data(), m.data(), m.size());
        }*/
        ManagedMemory(uint32_t s) { resize(s); }

        void Zero() const
        {
            std::memset((void*)data(), 0, size());
        }

        template < typename Y > Y & rt(uint32_t index=0, uint32_t g_index = 0)
        {
            return *pt<Y>(index,g_index);
        }

        template < typename Y > Y & rlt(uint32_t index=0)
        {
            return *plt<Y>(index);
        }

        template < typename T > T * pt(uint32_t index=0,uint32_t offset=0)
        {
            return (((T*)(data()+offset)) + index);
        }

        template < typename T > T * plt(uint32_t index=0)
        {
            return (((T*)(data()+size()))- (1+ index));
        }

        template <typename T, typename ... t_args> T * Allocate(t_args ... args)
        {
            uint32_t s = (uint32_t)size();
            resize(s+sizeof(T));

            new((data()+s)) T(args...);

            return (T*)(data()+s);
        }

        template <typename T, typename ... t_args> uint32_t AllocateOffset(t_args ... args)
        {
            uint32_t s = (uint32_t)size();
            resize(s+sizeof(T));

            new((data()+s)) T(args...);

            return s;
        }

        template <typename T > T* Offset(uint32_t o)
        {
            return (T*)(data()+o);
        }

        template <typename T>void append(const T &t)
        {
            auto os = size();
            resize(os + t.size());
            std::memcpy(data() + os, t.data(), t.size());
        }

        void Stream(uint32_t s)
        {
            resize(size() + sizeof(uint32_t));
            *plt<uint32_t>() = s;
        }

        /*void Stream(Memory m)
        {
            uint32_t offset = (uint32_t)size();
            resize(offset + m.size());
            std::memcpy(data()+offset,m.data(),m.size());
        }*/
    };

    /*
        A DerivativeIdentify takes two frequencies and generates a matrix solution which isolates each waveform.
    */

    template < typename T > class DerivativeIdentity
    {
        struct N
        {
            N(T * _p):p(_p){}

            T& operator[](T i) const 
            {
                return *(p+i);
            }

            T * p;
        };
    public:
        DerivativeIdentity(T _x, T _y)
            : x(_x)
            , y(_y)
        {
            if(x > y)
            {
                T t = x;
                x = y;
                y = t;
            }

            ManagedMemory tmp;
            tmp.resize(y * y * sizeof(T));
            T * p = tmp.pt<T>();

            auto abs = [](T v)-> T { if(v > 0) return v; return v*-1;};

            T d = y - x;
            //k = d + 1;
            T positive = x;
            T negative = 1;

            T pitr = positive;
            T nitr = negative;

            T itr = 2;
            p[0]=positive;
            p[1]=negative;
            while(d != abs(pitr-nitr))
            {
                pitr = (pitr + x) % y;

                nitr = nitr - x;
                if(nitr < 0)
                    nitr+=y;

                /*if(nitr == pitr)
                {
                    p[itr] = pitr;
                    itr++;
                    break;
                }

                if(nitr == positive || nitr == negative || pitr == negative || pitr == positive)
                    break;*/

                p[itr+1] = nitr;
                p[itr] = pitr;

                itr += 2;
            }

            sl = itr;

            m.resize(y*sl*sizeof(T));
            T *mp = m.pt<T>();

            for(T i=0;i<y;i++)
                for(T j = 0; j < sl; j++)
                    *mp++ = (p[j] + i) % y;
        }

        N operator[](T i) const
        {
            return N(((T*)m.data())+i*sl);
        }

    //private:
        //T d;
        //T k;
        T sl;
        T x;
        T y;
        ManagedMemory m;
    };

    /*
        A Rational Wave wraps an sample buffer and provides methods to comprehend the data.
    */

    template <typename T > class RationalWave
    {
    public:
        RationalWave(){}

        template <typename ... F>RationalWave(T rep, F ... args)
        {
            Generate(rep,args...);
        }

        template <uint64_t O> void _Generate(){}

        template <uint64_t O, typename Q, typename ... F> void _Generate(Q c, F ... args)
        {
            self_t[O] = (T)c;
            _Generate<O+1>(args...);
        }

        template <typename ... F> void Generate(T rep, F ... args)
        {
            constexpr T unit =  sizeof...(args);
            set_size(unit*rep);

            _Generate<0>(args...);

            T gap = unit * sizeof(T);
            for(T i = 1; i < rep; i++)
                std::memcpy(m.data() + gap * i,m.data() ,gap);
        }

        /*void Generate2(T r, T s1, T s2)
        {
            set_size(2*r);
            T * p = m.pt<T>();

            for(T i = 0;i<r;i++)
            {
                p[i*2] = s1;
                p[i*2 + 1] = s2;
            }
        }

        void Generate3(T r, T s1, T s2, T s3)
        {
            set_size(3*r);
            T * p = m.pt<T>();

            for(T i = 0;i<r;i++)
            {
                p[i*3] = s1;
                p[i*3 + 1] = s2;
                p[i*3 + 2] = s3;
            }
        }*/

        template < typename F > void PairIterator(RationalWave & r, F f)
        {
            if(r.m.size() != m.size())
                throw 1;

            T * p1 = m.pt<T>();
            T * p2 = r.m.pt<T>();

            T l = (T)(m.size()/sizeof(T));
            while(l--)
                f(*p1++,*p2++);		
        }

        void operator -=(RationalWave & r)
        {
            PairIterator(r,[](T & l, T & r) { l -= r; });
        }

        void operator +=(RationalWave & r)
        {
            PairIterator(r,[](T & l, T & r) { l += r; });	
        }

        T & operator[](T i){ return *m.pt<T>((uint32_t)i); }

        T samples() { return (T)(m.size()/sizeof(T)); }
        T * root() { return m . pt<T>(); }
        void set_size(T c)
        {
            m.resize(c*sizeof(T));
        }

        RationalWave DerivativeFrequency(T frequency)
        {
            RationalWave rw;

            T l = (T)((m.size()/sizeof(T))/frequency);
            if((m.size()/sizeof(T))%frequency)
                throw 1;

            T sz = samples();
            rw.set_size(sz);

            for(T i = 0; i<sz;i++)
                rw[(i + l) % sz] = self_t[(i + l) % sz] - self_t[i];

            return rw;
        }

        void Multiply(T m)
        {
            auto r = root();
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

            T sz = samples();
            if(sz%d.y)
                throw 1;

            T r = sz/d.y;

            rw.set_size(sz);

            for(T i = 0; i < d.y; i++)
            {
                for(T j = 0; j<d.sl;j++)
                {
                    if(j==0)
                        rw[i] = self_t[d[i][j]];
                    else
                        rw[i] += self_t[d[i][j]];
                }
            }

            for(T i = 1; i < r; i++ )
                std::memcpy(rw.m.data() + i *d.y*sizeof(T) ,rw.m.data(), d.y * sizeof(T));

            return rw;
        }

        RationalWave AntiDerivative(T amplitude)
        {
            RationalWave rw;

            T sz = samples();
            rw.set_size(sz);

            rw[0] = amplitude;

            for(T i = 1; i < sz; i++)
            {
                amplitude += self_t[i-1];
                rw[i] = amplitude;
            }

            return rw;
        }

        std::pair<bool,bool> IsBalancedAndRegular(T frequency)
        {
            bool balance = true;
            bool regular = true;
            T l = (T)((m.size()/sizeof(T))/frequency);
            if((m.size()/sizeof(T))%frequency)
                throw 1;

            T sz = samples();
            T sum = 0;

            for(T i = 0; i < sz; i++)
            {
                if(self_t[i] != self_t[(i + frequency)%sz])
                    regular = false;

                if(i % frequency == 0)
                {
                    if(sum)
                    {
                        balance = false;
                        sum = 0;
                    }
                }
                sum += self_t[i];
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
            auto r = root();
            
            T l = samples();

            for(T i = 0; i < l;i++)
                std::cout << (int64_t)r[i] << " ";
            
            std::cout << std::endl << std::endl;
        }

    private:
        ManagedMemory m;
    };
}