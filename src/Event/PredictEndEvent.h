#ifndef SOFA_SIMULATION_PREDICTENDEVENT_H
#define SOFA_SIMULATION_PREDICTENDEVENT_H

#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/config.h>
#include <string>

class SOFA_SIMULATION_CORE_API PredictEndEvent : public sofa::core::objectmodel::Event
{
public:

virtual size_t getEventTypeIndex() const override { return PredictEndEvent::s_eventTypeIndex; }
static bool checkEventType( const Event* event ) { return event->getEventTypeIndex() == PredictEndEvent::s_eventTypeIndex; }
virtual const char* getClassName() const override {return name.c_str();}

PredictEndEvent( SReal dt);
PredictEndEvent() {}
void setDt(SReal sdt){dt = sdt;}
~PredictEndEvent() override;

SReal getDt() const { return dt; }
inline static const char* GetClassName() { return "PredictEndEvent"; }

protected:
SReal dt;
static const size_t s_eventTypeIndex;
std::string name;
};
#endif
