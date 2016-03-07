#include "Debugger.h"
Debugger::Debugger(){
  this->m_debug = false;
  this->pfout = NULL;
}
Debugger::Debugger(std::string filename, bool set){
  this->setFileName(filename);
  this->setDebug(set);
}
Debugger::~Debugger(){
  this->clear();
}
void Debugger::clear(){
  if(this->pfout != NULL){
    if(this->pfout->is_open()){
      this->pfout->close();
    }
    delete pfout;
    pfout = NULL;
  }
}
void Debugger::create(){
  this->pfout = new std::ofstream();
}
bool Debugger::IsSet(){
  return this->m_debug;
}
void Debugger::setDebug(bool debug){
  this->m_debug = debug;
  if(debug == true){
    this->create();
  }else{
    this->clear();
  }
}
void Debugger::setFileName(const std::string& filename){
 this->filename = filename;
}
void Debugger::debug(const std::string& message){
  if(this->m_debug){
      this->pfout->open(this->filename.c_str(), std::ofstream::app);
      *(this->pfout) << message << std::endl;
      this->pfout->close();
    }
}
