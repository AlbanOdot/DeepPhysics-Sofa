#include "PredictEndEvent.h"



SOFA_EVENT_CPP( PredictEndEvent )

PredictEndEvent::PredictEndEvent(SReal dt, int id=-1)
: sofa::core::objectmodel::Event()
, dt(dt)
{
    if(id < 0)
        name = std::string("PredictEndEvent");
    else
        name = std::string("PredictEndEvent_").append(std::to_string(id)).c_str();
}


PredictEndEvent::~PredictEndEvent()
{
}
