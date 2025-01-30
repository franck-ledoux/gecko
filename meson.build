project('gecko', 'cpp',
        version: '1.0',
        default_options: ['warning_level=3', 'cpp_std=c++17'])

#==========================================================
# dependencies to extern projects
#==========================================================
boost = subproject('boost-1.87.0')
boost_dep = boost.get_variable('boost_dep')

cgal = subproject('CGAL-5.6.2')
cgal_dep = cgal.get_variable('cgal_dep')

predicates = subproject('predicates')
predicates_dep = predicates.get_variable('predicates_lib_dep')

eigen3_dep = dependency('eigen3')

catch2_proj = subproject('catch2')
catch2_dep = catch2_proj.get_variable('catch2_dep')
catch2_main_dep = catch2_proj.get_variable('catch2_with_main_dep')

nlohmann_json_dep = dependency('nlohmann_json')
#==========================================================
# path to the test repo
#==========================================================
test_samples_dir = '"'+meson.project_source_root()+'/test_data"'
configure_file(
        output: 'unit_test_config.h',
        configuration: configuration_data({ 'TEST_SAMPLES_DIR': test_samples_dir}),
        install: false
)
# denpendency tyo provide to unit test suites
ut_dep = declare_dependency(include_directories : '.')

#==========================================================
# Project structure
#==========================================================
subdir('gmds_core')
subdir('gblock')
subdir('mcts')
