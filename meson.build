project('gecko', 'cpp',
        version: '1.0',
        default_options: ['warning_level=3', 'cpp_std=c++14'])



cgal_dep = dependency('cgal', required : false)
eigen3_dep = dependency('eigen3')

catch2_proj = subproject('catch2')
catch2_dep = catch2_proj.get_variable('catch2_dep')
catch2_main_dep = catch2_proj.get_variable('catch2_with_main_dep')


# path to the test repo
test_samples_dir = 'path/to/test/samples'
configure_file(
        output: 'unit_test_config.h',
        configuration: configuration_data({ 'TEST_SAMPLES_DIR': test_samples_dir}),
        install: false
)


subdir('externals')
subdir('gmds_core')
subdir('gblock')
subdir('tst')
