#include "PredictionPickedEvent.h"

SOFA_EVENT_CPP( PredictionPickedEvent )

PredictionPickedEvent::PredictionPickedEvent(SReal dt, int id=-1)
: sofa::core::objectmodel::Event()
, dt(dt)
{
    if(id < 0)
        name = std::string("PredictionPickedEvent");
    else
        name = std::string("PredictionPickedEvent_").append(std::to_string(id)).c_str();
}


PredictionPickedEvent::~PredictionPickedEvent()
{
}
