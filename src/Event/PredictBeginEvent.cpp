#include "PredictBeginEvent.h"



SOFA_EVENT_CPP( PredictBeginEvent )

PredictBeginEvent::PredictBeginEvent(SReal dt)
        : sofa::core::objectmodel::Event()
        , dt(dt)
{
        name = std::string("PredictBeginEvent");
}


PredictBeginEvent::~PredictBeginEvent()
{
}

