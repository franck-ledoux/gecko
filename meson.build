project('gecko', 'cpp',
        version: '1.0',
        default_options: ['warning_level=3', 'cpp_std=c++17'])



cgal_dep = dependency('cgal', required : false)
eigen3_dep = dependency('eigen3')

catch2_proj = subproject('catch2')
catch2_dep = catch2_proj.get_variable('catch2_dep')
catch2_main_dep = catch2_proj.get_variable('catch2_with_main_dep')


# Sous-répertoires
subdir('externals')
subdir('gmds_core')
subdir('gblock')
subdir('tst')
