
set(VERSION "1.0.0")
add_definitions("-DVERSION=\"${VERSION}\"")

include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(NO_X_SRCS oshook.cxx sequoia.cxx)
set(WITH_X_SRCS xoshook.cxx xsequoia.cxx mtwin.c seqwin.c matwin.c strwin.c cvector.c)

set(MATRIX_SRCS array.cxx vector.cxx matrix.cxx rmat.cxx sequence.cxx mystring.cxx
        brookhav.cxx path_mat.cxx parse.cxx overlay.cxx comp_mat.cxx ltable.cxx parse2.cxx
        instanti.cxx instant2.cxx)

#object files in common
set(BOTH_SRCS ${MATRIX_SRCS} consens.cxx readline.c
	geoscore.cxx timelice.cxx sort.cc
	variable.cxx)

set(CK_SRCS ${BOTH_SRCS} ${NO_X_SRCS})

set(CK_HDRS array.hxx array2d.hxx matrix.hxx sequence.hxx 
	comp_mat.hxx path_mat.hxx parse.hxx mystring.hxx 
	vector.hxx brookhav.hxx overlay.hxx rmat.hxx 
	rvec.hxx cmbmacro.hxx consens.hxx oshook.hxx variable.hxx)

set(XCK_SRCS ${BOTH_SRCS} ${WITH_X_SRCS})

set(XCK_HDRS array.hxx array2d.hxx matrix.hxx sequence.hxx 
	comp_mat.hxx path_mat.hxx parse.hxx mystring.hxx 
	vector.hxx brookhav.hxx overlay.hxx rmat.hxx 
	rvec.hxx oshook.hxx cmbmacro.hxx consens.hxx xsequoia.hxx 
	mtwin.h)

########################
### READLINE LIBRARY ###
########################

find_library(READLINE_LIBRARY libreadline.a
    PATHS "C:/Program Files/GnuWin32/lib")
set(LIBS ${LIBS} ${READLINE_LIBRARY})
if(MSVC)
    find_library(GNUWIN32_LIBRARY libgw32c.a
        PATHS "C:/Program Files/GnuWin32/lib")
    set(LIBS ${LIBS} ${GNUWIN32_LIBRARY})
endif()
find_path(READLINE_INCLUDE_DIR readline/readline.h
    PATHS "C:/Program Files/GnuWin32/include")
if(READLINE_INCLUDE_DIR)
    include_directories(${READLINE_INCLUDE_DIR})
endif()
add_definitions("-DREADLINE_STATIC")
find_library(NCURSES_LIBRARY ncurses)
if(NCURSES_LIBRARY)
    set(LIBS ${LIBS} ${NCURSES_LIBRARY})
endif()


add_executable(sequoia ${CK_SRCS} ${CK_HDRS})
target_link_libraries(sequoia ${LIBS})

set(SEQUOIA_BUILD_XSEQUOIA OFF CACHE BOOL "Whether to build X-windows version of sequoia")
if(SEQUOIA_BUILD_XSEQUOIA)
    add_executable(xsequoia ${XCK_SRCS} ${XCK_HDRS})
    target_link_libraries(xsequoia ${LIBS})
endif()

site_name(HOSTNAME)
if(WIN32)
    execute_process( COMMAND "date /T"
        OUTPUT_VARIABLE DATETIME 
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_VARIABLE DATETIME 
        ERROR_STRIP_TRAILING_WHITESPACE
    )
    set(BUILD_USER "$ENV{USERNAME}")
else()
    execute_process( COMMAND date OUTPUT_VARIABLE DATETIME OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process( COMMAND whoami OUTPUT_VARIABLE BUILD_USER OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

configure_file(timestamp.h.in ${CMAKE_CURRENT_BINARY_DIR}/timestmp.h @ONLY)

add_executable(vartest vartest.cxx ${BOTH_SRCS} oshook.cxx)
target_link_libraries(vartest ${LIBS})
add_test(vartest vartest)

