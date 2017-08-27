//          Copyright Federico Kircheis 2017
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef CALCNUM_BISEC_HPP_0925487529
#define CALCNUM_BISEC_HPP_0925487529

#include "intervals.hpp"
#include "utils.hpp"

#include <functional>
#include <cassert>
namespace calcnum{

    // length(interval) and max(|fa|,|fb|), min(|fa|,|fb|) can be used as a stopping criteria
    struct bisect_status{
        closed_interval interval;
        double fa;
        double fb;
    };

    class bisect_iter{
        std::function<double(double)> f;
        bisect_status st;

    public:
        bool invariant() const {
            return (signum(st.fa)*signum(st.fb) <= 0) && calcnum::invariant(st.interval);
        }

        explicit bisect_iter(const std::function<double(double)>& f_, const closed_interval& interval) : f(f_){
            st.interval = interval;
            st.fa = f(interval.a);
            st.fb = f(interval.b);

            test_invariant<bisect_iter> _(*this);
        }

        bisect_iter& operator++(){
            test_invariant<bisect_iter> _(*this);
            const auto mid = midpoint(st.interval);
            const auto fmid = f(mid);

            const auto signfmid = signum(fmid);
            const auto signfa = signum(st.fa);
            const auto signfb = signum(st.fb);
            if(signfmid == 0){ // just in case we have really hit the 0, will never happen with real data...
                st.interval = {mid, mid};
                st.fa = fmid;
                st.fb = fmid;
            } else if(signfmid*signfa == sign::greaterzero){
                st.interval.a = mid;
                st.fa = fmid;
            } else if(signfmid*signfb == sign::greaterzero){
                st.interval.b = mid;
                st.fb = fmid;
            } else {
                assert(false && "FIXME: assert if some element is nan, otherwise should never happen");
            }
            return *this;
        }

        const bisect_status& operator*() const {
            test_invariant<bisect_iter> _(*this);
            return st;
        }
        const bisect_status* operator->() const {
            test_invariant<bisect_iter> _(*this);
            return &st;
        }
    };
}

#endif
