#include "my_class.h"
#include <iostream>

MyClass::MyClass(const std::string& name) : name_(name) {}

void MyClass::say_hello() const {
    std::cout << "Hello, " << name_ << "!" << std::endl;
}
