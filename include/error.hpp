#ifndef CALCNUM_ERROR_HPP_8826792364
#define CALCNUM_ERROR_HPP_8826792364

#include <cfenv>
#include <cmath>
#include <cerrno>
#include <cassert>


namespace calcnum{

#if (math_errhandling & MATH_ERRNO)
	struct reset_errno {
		const int m_errno; // verify, should be int, but compiler may annotate it (should therefore prefer decltype?)
		reset_errno() noexcept: m_errno(errno)  {
			errno = {};
		}
		~reset_errno(){
			errno = m_errno;
		}
	};
#endif

	inline std::fenv_t get_fenv(){
		std::fenv_t envp{};
		auto res = std::fegetenv(&envp);
		assert(res==0 && "FIXME: handle possible failure of fegetenv");
		return envp;
	}

#if (math_errhandling & MATH_ERREXCEPT)
    struct reset_fenv {
        const std::fenv_t m_envp;
        reset_fenv() : m_envp(get_fenv()){
            auto res = std::feclearexcept(FE_ALL_EXCEPT);
            assert(res==0 && "FIXME: handle possible failure of feclearexcept");
        }
        ~reset_fenv() {
            const auto res = std::fesetenv(&m_envp);
            assert(res==0 && "FIXME: handle possible failure of fesetenv");
        }
    };
#endif



    struct reset_float_env {
#if (math_errhandling & MATH_ERRNO)
        reset_errno err;
#endif
#if (math_errhandling & MATH_ERREXCEPT)
		reset_fenv fenv;
#endif
	};

}
#endif
