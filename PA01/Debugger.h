#ifndef __DEBUGGER__H__
#define __DEBUGGER__H__
#include <fstream>
class Debugger{
  protected:
    bool m_debug;
    std::string filename;
    std::ofstream *pfout;

  public:
    Debugger();
    Debugger(std::string filename, bool set);
    ~Debugger();
    void create();
    void clear();
    bool IsSet();
    void setDebug(bool debug);
    void setFileName(const std::string& message);
    void debug(const std::string& message);
};
#endif