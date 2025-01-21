project('gecko', 'cpp',
        version: '1.0',
        default_options: ['cpp_std=c++17'])

# Dépendances externes éventuelles (exemple avec Boost)
boost_dep = dependency('boost', required: false) # Facultatif

# Sous-répertoires
subdir('lib')
subdir('src')


executable('gecko',
           files('main.cpp'),
           include_directories: include_directories('../include'),
           link_with: gecko_lib,  # Lier avec la bibliothèque libmy_library.a
)
