#ifndef APP_H
#define APP_H

#include <string>

struct Options {
    std::string path;
    int loadNframes;
    std::string encoding;
    bool do_normalize;
};

class App {
private:
    Options* options;
public:
    void run();
    App(Options &options);
    ~App();
};

#endif // APP_H
