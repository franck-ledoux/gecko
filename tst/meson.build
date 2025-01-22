tst_files = files('test2.cpp')

catch2_proj = subproject('catch2')
catch2_dep = catch2_proj.get_variable('catch2_dep')
catch2_main_dep = catch2_proj.get_variable('catch2_with_main_dep')



tst_exe = executable('gecko2_test',
                      tst_files,
                     dependencies : [ catch2_dep, catch2_main_dep ])

test('gecko2_test', tst_exe)