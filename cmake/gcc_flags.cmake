##################################################################
# Settings for the GCC compiler (should work with >=4.8)

##################################################################
# language, use CMAKE_CXX_STANDARD, CMAKE_CXX_STANDARD_REQUIRED, CXX_EXTENSIONS if possible
#add_compile_options(-std=c++14)

##################################################################
# compiler warnings, should be enabled on every project
# unlike msvc, -Wall does not enable all warnings

# generic warnings/errors
add_compile_options(
	-Wall
	-Wextra
	-pedantic-errors
	-Werror=pedantic
	-Werror=main
	-Wunreachable-code
	-Wunused
	-Wunknown-pragmas
	-Wpragmas
	-Werror=return-type
)

# portability
add_compile_options(
	-Werror=multichar
	-Werror=address
	-Werror=sequence-point
	-Werror=cpp
	-Werror=strict-aliasing
	-Werror=strict-null-sentinel
	-Werror=trigraphs
)

# extensions
add_compile_options(
	-Werror=pointer-arith
	-fno-nonansi-builtins
	-fstrict-enums
	-fvisibility-inlines-hidden
)

# multiple declaration, shadowing, eval undefined identifier
add_compile_options(
	-Wshadow
	-Wundef
	-Werror=redundant-decls
)

# others
add_compile_options(
	-Waggressive-loop-optimizations # generic UB
	-Werror=free-nonheap-object     # ponter arith
	-Wignored-qualifiers
	-Wmissing-declarations   # should set function in anonym namespace or declare it to avoid possible clashes
	-Wnarrowing
	-Wpacked
	-Wparentheses
	-Werror=return-local-addr

	# unused /unnecessary code
	-Wunused-but-set-parameter
	-Wunused-but-set-variable
	-Wunused-function
	-Wunused-label
	-Wunused-macros
	-Wunused-parameter
	-Wunused-result
	-Wunused-value
	-Wunused-variable
	-Wuseless-cast

	-Wvarargs
)

if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 5)
	add_compile_options(
		-Wodr
		-Wsized-deallocation
	)
	if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 6)
		add_compile_options(
			-Wignored-attributes
			-Wmisleading-indentation
			-Wsubobject-linkage
			-Wunused-const-variable=2
		)
		if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 7)
			add_compile_options(
				-Wbuiltin-declaration-mismatch
			)
		endif()
	endif()
endif()

# macro
add_compile_options(
	-Werror=builtin-macro-redefined
	-Werror=endif-labels
)
if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 7)
	add_compile_options(
		-Wexpansion-to-defined
	)
endif()

# deprecation
add_compile_options(
	-Wdeprecated
	-Wdeprecated-declarations
)

# class/structs and init
add_compile_options(
	-Wctor-dtor-privacy
	-Werror=non-virtual-dtor
	-Werror=reorder
	-Werror=uninitialized
	-Werror=maybe-uninitialized
	-Werror=delete-non-virtual-dtor
	-Winherited-variadic-ctor
	-Werror=init-self
	-Werror=invalid-offsetof
	-Wmissing-field-initializers
	-Wvirtual-move-assign
)
if(CMAKE_CXX_COMPILER_VERSION LESS_EQUAL 4)
	add_compile_options( -Wno-missing-field-initializers) # otherwise S s = {}; issues warning
endif()

# switch/branches
add_compile_options(
	-Wswitch
	-Wswitch-default
	-Wswitch-enum
	-Wempty-body
	-Wenum-compare
)
if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 6)
	add_compile_options(
		-Wduplicated-cond
	)
	if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 7)
		add_compile_options(
			-Wdangling-else
			-Wduplicated-branches
			-Wimplicit-fallthrough
			-Wswitch-unreachable
		)
	endif()
endif()

# nullptr and casts warnings (may generate a lot of warnings when interacting with old code or C code)
add_compile_options(
	-Wzero-as-null-pointer-constant
	-Wold-style-cast
	-Wuseless-cast
	-Wnonnull
	-Wcast-qual
	-Wcast-align
	-Werror=null-dereference
	-Wconversion-null
	-Wint-to-pointer-cast
)


# arithmetic/numeric warnings
add_compile_options(
	-Wconversion
	-Wfloat-equal
	-Wsign-compare
	-Wsign-promo
	-Wsign-conversion
	-Werror=div-by-zero
	-ftrapv
)
if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 5)
	add_compile_options(
		-Werror=shift-count-negative
		-Werror=shift-count-overflow
	)
	if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 6)
		add_compile_options(
			-Wshift-negative-value
			-Werror=shift-overflow
		)
	endif()
endif()



# logical operations
add_compile_options(
	-Wlogical-op
)

# possible code structure problem
add_compile_options(
	-Wdisabled-optimization
)

# formatting
add_compile_options(
	-Wformat=2
	-Wformat-nonliteral
	-Wformat-y2k
	-Werror=format
	-Werror=format-security
	-Werror=format-extra-args
	-Wformat-contains-nul
	-Wformat-zero-length
)
if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 5)
	add_compile_options(
		-Wformat-signedness
	)
	if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 7)
		add_compile_options(
			-Wformat-overflow
			-Wformat-truncation
		)
	endif()
endif()

# exception safety
add_compile_options(
	-Werror=terminate
)

# operations on booleans
if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 5)
	add_compile_options(
		-Werror=bool-compare
		-Wlogical-not-parentheses
		-Wswitch-bool
	)

	if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 5)
		add_compile_options(
			-Wtautological-compare
		)
		if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 7)
			add_compile_options(
				-Werror=bool-operation
				-Wint-in-bool-context
			)
		endif()
	endif()
endif()
# suggestions for improving code
add_compile_options(
	# formatting
	-Wsuggest-attribute=format
	# function signature
	-Wsuggest-attribute=const
	-Wsuggest-attribute=pure
	-Wsuggest-attribute=noreturn
)
if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 5)
	add_compile_options(
		# inheritance
		-Wsuggest-final-methods
		-Wsuggest-final-types
		-Wsuggest-override
	)
endif()


# memory
add_compile_options(
	-Waddress
	-Werror=vla
	-Werror=array-bounds
	-Werror=write-strings
	-Werror=overflow
	-Wsizeof-pointer-memaccess
	-Wsizeof-array-argument
)
if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 5)
	add_compile_options(
		-Werror=memset-transposed-args
	)
	if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 6)
		add_compile_options(
			-Wchkp
			-Wplacement-new=2
			-Wscalar-storage-order
			-Wnull-dereference
		)
		if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 7)
			add_compile_options(
				-Waligned-new=all
				-Werror=alloc-zero
				-Werror=alloca
				-Werror=stringop-overflow
				-Wmemset-elt-size
				-Wpointer-compare
			)
		endif()
	endif()
endif()

option(EFFCXX "Enable Effective C++ Warnings (very noisy on correct c++11 code)"  OFF)
if(EFFCXX)
	message(" ===== Enabled Effective C++ Warnings =====")
	add_compile_options(-Weffc++)
endif()


##################################################################
# project structure
add_compile_options(
	-Wmissing-include-dirs
)

##################################################################
# additional debug informations
option(DEBUG_ITERATORS "Additional debug information with glibc" OFF)
if(DEBUG_ITERATORS)
	message(" ===== Enabled additional debug information for iterators =====")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC")
endif()


##################################################################
# sanitizers (checks made at runtime)
# You can choose only one(!)
set(SanitizerValues "NONE;SANADD;SANUNDEF;SANTHREAD" CACHE STRING
	"List of possible values for the sanitizers")
set(SanValue "SANUNDEF" CACHE STRING
	"SanValue chosen by the user at CMake configure time")
set_property(CACHE SanValue PROPERTY STRINGS ${SanitizerValues})

# not using add_compile_options since since flag is also used by the linker
if("${SanValue}" STREQUAL "SANADD")
	message(" ===== Enable sanitizer for addresses ===== ")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
elseif("${SanValue}" STREQUAL "SANUNDEF")
	message(" ===== Enable sanitizer for UB ===== ")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined")
elseif("${SanValue}" STREQUAL "SANTHREAD")
	message(" ===== Enable sanitizer for threads ===== ")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread")
endif()

option(PROFILE "Enable profiling" OFF)
if(PROFILE)
	message(" ===== Enable profiling ===== ")
	add_compile_options(-pg)
endif()

##################################################################
# mitigations
# abort and provide some useful message instead of ignoring errors, crash at another random place or leave some vulnerability open
option(FORTIFY "Enable fortify sources (already enabled for release target, also enables optimizations)" OFF)
if(FORTIFY)
	message(" ===== Enable fortify sources (always enabled for release target, also enables optimizations) ===== ")
	add_definitions(-D_FORTIFY_SOURCE=2 -O2)
endif()
# FIXME: add check if compiling with -O1, in that case we should use -D_FORTIFY_SOURCE=1...
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -D_FORTIFY_SOURCE=2") # =1 when using -O1..

# not using add_compile_options since since flag is also used by the linker
option(RELRO "Enable full relro" OFF)
if(RELRO)
    message(" ===== Enabled full relro =====")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,-z,relro,-z,now")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,-z,relro")
endif()

option(STACK_PROTECTOR "Enable stack protector on all functions")
if(STACK_PROTECTOR)
    message(" ===== Enabling full stack protector =====")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fstack-protector-all")
else()
    if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 5)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fstack-protector-strong")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fstack-protector")
    endif()
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fstack-check")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pie -fpie -fpic -fPIC")

if( CMAKE_SYSTEM_NAME STREQUAL "Windows" )
	# enable ASLR
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,dynamicbase ")
	# enable DEP
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,nxcompat")
endif()

option(CHANGE_CHAR_SIGN "Change sign of char to test code for portability"  OFF)
if(CHANGE_CHAR_SIGN)
	set(sign_char_values "-funsigned-char;-fsigned-char" CACHE STRING "List of possible values for the sign of char cache variable")
	set(sign_char "-fsigned-char" CACHE STRING "Value chosen by the user at configure time")
	set_property(CACHE sign_char PROPERTY STRINGS ${sign_char_values})
	add_compile_options(${sign_char})
endif()

