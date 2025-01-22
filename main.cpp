#include <iostream>
#include "my_class.h"

int main() {
    std::cout << "Hello, Meson + C++!" << std::endl;

    MyClass obj("World");
    obj.say_hello();

    return 0;
}
