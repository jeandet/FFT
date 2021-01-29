project('FFT', 'cpp',default_options : ['cpp_std=c++2a'])


FFT_inc = include_directories(['include'])

FFT_lib = library('FFT' ,extra_files:'include/FFT/FFT.hpp')

FFT_dep = declare_dependency(include_directories: FFT_inc)

catch2_dep = dependency('catch2', main : true, fallback : ['catch2', 'catch2_dep'])

foreach test:['simple']
    exe = executable(test,'tests/'+test+'/main.cpp',
                    dependencies:[catch2_dep, FFT_dep],
                    install: false
                    )
    test(test, exe)
endforeach
