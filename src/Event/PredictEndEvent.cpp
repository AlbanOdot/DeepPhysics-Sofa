#include "PredictEndEvent.h"



SOFA_EVENT_CPP( PredictEndEvent )

PredictEndEvent::PredictEndEvent(SReal dt)
: sofa::core::objectmodel::Event()
, dt(dt)
{
    name = std::string("PredictEndEvent");
}


PredictEndEvent::~PredictEndEvent()
{
}
