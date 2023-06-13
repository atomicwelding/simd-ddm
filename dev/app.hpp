#ifndef APP_H
#define APP_H

#include <string>

class App {
private:
    const std::string& path;
public:
    void run();
    App(const std::string& path);
    ~App();
};

#endif // APP_H
