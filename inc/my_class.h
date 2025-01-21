#ifndef MY_CLASS_H
#define MY_CLASS_H

#include <string>

class MyClass {
private:
    std::string name_;

public:
    MyClass(const std::string& name);
    void say_hello() const;
};

#endif // MY_CLASS_H
