#ifndef SOFA_SIMULATION_PREDICTBEGINEVENT_H
#define SOFA_SIMULATION_PREDICTBEGINEVENT_H

#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/config.h>
#include <string>

class SOFA_SIMULATION_CORE_API PredictBeginEvent : public sofa::core::objectmodel::Event
{
public:


virtual size_t getEventTypeIndex() const override { return PredictBeginEvent::s_eventTypeIndex; }
static bool checkEventType( const Event* event ) { return event->getEventTypeIndex() == PredictBeginEvent::s_eventTypeIndex; }
virtual const char* getClassName() const override {return name.c_str();}

PredictBeginEvent( SReal dt);
PredictBeginEvent() {}
void setDt(SReal sdt){dt = sdt;}
~PredictBeginEvent() override;

SReal getDt() const { return dt; }
inline static const char* GetClassName() { return "PredictBeginEvent";}
protected:
SReal dt;
static const size_t s_eventTypeIndex;
std::string name;
};

#endif
