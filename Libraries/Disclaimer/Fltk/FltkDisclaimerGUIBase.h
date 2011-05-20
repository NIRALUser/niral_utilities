#ifndef FltkDisclaimerGUIBase_H
#define FltkDisclaimerGUIBase_H

class FltkDisclaimerGUIBase
{
public:
  virtual ~FltkDisclaimerGUIBase() {};
  virtual void Accept() = 0;
  virtual void Reject() = 0;
};
  
#endif
