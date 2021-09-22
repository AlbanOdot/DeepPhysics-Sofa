#pragma once
#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/config.h>
#include <string>

class SOFA_SIMULATION_CORE_API PredictionPickedEvent : public sofa::core::objectmodel::Event
{
public:

virtual size_t getEventTypeIndex() const override { return PredictionPickedEvent::s_eventTypeIndex; }
static bool checkEventType( const Event* event ) { return event->getEventTypeIndex() == PredictionPickedEvent::s_eventTypeIndex; }
virtual const char* getClassName() const override {return name.c_str();}

PredictionPickedEvent( SReal dt, int id );
PredictionPickedEvent() {}
void setDt(SReal sdt){dt = sdt;}
~PredictionPickedEvent() override;

SReal getDt() const { return dt; }
inline static const char* GetClassName() { return "PredictionPickedEvent"; }

inline  uint getID() { return this->id;}
inline void setID(const uint id) {this->id = id;}
protected:
SReal dt;
uint id;
static const size_t s_eventTypeIndex;
std::string name;
};
