#include "PredictBeginEvent.h"



SOFA_EVENT_CPP( PredictBeginEvent )

PredictBeginEvent::PredictBeginEvent(SReal dt, int id=-1)
        : sofa::core::objectmodel::Event()
        , dt(dt)
{
    if(id < 0)
        name = std::string("PredictBeginEvent");
    else
        name = std::string("PredictBeginEvent_").append(std::to_string(id)).c_str();
}


PredictBeginEvent::~PredictBeginEvent()
{
}

