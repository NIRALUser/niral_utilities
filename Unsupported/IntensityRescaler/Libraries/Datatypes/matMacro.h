#ifndef _matMacro_H 
#define _matMacro_H 

#define matSetMacro(name,type) \
  type m_##name; \
    virtual void Set##name (type _arg) \
{ \
    this->m_##name = _arg; \
} 

#define matGetMacro(name,type) \
    virtual type Get##name () const { \
    return this->m_##name; \
} 

#define matSetFilterMacro(caption,name,type,defaultparam) \
  m_param.param = #name ; \
  m_param.paramcaption = #caption; \
  m_param.paramtype = #type; \
  type* m_##name = new type(); \
  m_param.pointer = (void*)m_##name;  \
  m_param.paramdefault = #defaultparam; \
  m_list.push_back(m_param);


#define matApplyFilterMacro(function,name,value) \
  function ->Set##name ( value );

#define matSetupMacro(name,type) \
 if (m_name == #name)\
    filter->Set##name ( *((type *)m_value) );


#define ABS(x) (((x) > 0) ? (x) : (-1*(x)))

#endif 
